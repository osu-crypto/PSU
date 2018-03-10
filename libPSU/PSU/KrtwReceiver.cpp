#include "KrtwReceiver.h"

#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Crypto/Commit.h>
#include <cryptoTools/Network/Channel.h>
#include "libPSU/PsuDefines.h"
#include "Tools/SimpleIndex.h"

using namespace std;

namespace osuCrypto
{
	void KrtwReceiver::init(u64 psiSecParam, PRNG & prng, span<block> inputs, span<Channel> chls)
	{
		mPsiSecParam = psiSecParam;
		mPrng.SetSeed(prng.get<block>());

		senderOprf.configure(false, psiSecParam, 128);
		
		u64 baseCount= senderOprf.getBaseOTCount();
		mBaseChoice.resize(baseCount);
		mBaseChoice.randomize(mPrng);
		mBaseOTRecv.resize(baseCount);
		NaorPinkas baseOTs;
		baseOTs.receive(mBaseChoice, mBaseOTRecv, mPrng, chls[0], 1);
		senderOprf.setBaseOts(mBaseOTRecv, mBaseChoice);

	}
	void KrtwReceiver::output(span<block> inputs, span<Channel> chls)
	{

			u64 numThreads(chls.size());

		SimpleIndex simple;
		simple.init(inputs.size(),true);
		simple.insertItems(inputs, numThreads);
		//simple.print();

		std::cout << IoStream::lock << "Receiver: " << simple.mMaxBinSize << "\t " <<simple.mNumBins
			<< std::endl << IoStream::unlock;


		u64 theirMaxBinSize = simple.mMaxBinSize - 1; //assume same set size, sender has mMaxBinSize, receiver has mMaxBinSize+1
		u64	numOTs = simple.mNumBins*(theirMaxBinSize);

		std::vector<block> coeffs;
		
		std::vector<std::vector<block>> Ss(simple.mNumBins);
		for (u64 i = 0; i < simple.mNumBins; i++)
		{
			Ss[i].resize(theirMaxBinSize);
			for (u64 j = 0; j < theirMaxBinSize; j++)
				Ss[i][j] = mPrng.get<block>();
		}

		std::cout << IoStream::lock << "mBins[1].items[1]  " << simple.mBins[1].items[1] << std::endl << IoStream::unlock;
		std::cout << IoStream::lock << "Ss[1] " << Ss[1][1]<< std::endl << IoStream::unlock;

		senderOprf.init(2*numOTs, mPrng, chls[0]);

		std::cout << IoStream::lock << "senderOprf.init done" << std::endl << IoStream::unlock;

		//poly
		u64 polyMaskBytes = 128 / 8; //(mPsiSecParam + log2(simple.mMaxBinSize + 1) + 7) / 8;

		std::cout << IoStream::lock << "polyMaskBytes " << polyMaskBytes << std::endl << IoStream::unlock;


		auto routine = [&](u64 t)
		{
			auto& chl = chls[t];
			u64 binStartIdx = simple.mNumBins * t / numThreads;
			u64 tempBinEndIdx = (simple.mNumBins * (t + 1) / numThreads);
			u64 binEndIdx = std::min(tempBinEndIdx, simple.mNumBins);
			
			std::cout << IoStream::lock;
			polyNTL poly;
			poly.NtlPolyInit(polyMaskBytes);//length=lambda +log(|Y|)
			std::cout << IoStream::unlock;

			for (u64 i = binStartIdx; i < binEndIdx; i += stepSize)
			{
				auto curStepSize = std::min(stepSize, binEndIdx - i);

				senderOprf.recvCorrection(chl, curStepSize*theirMaxBinSize);

				std::vector<u8> sendBuff(curStepSize*theirMaxBinSize*(simple.mMaxBinSize + 1)*polyMaskBytes);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{

						std::vector<block> setY(simple.mBins[binIdx].mBinRealSizes, Ss[binIdx][itemTheirIdx]);
						std::vector<block>sendEncoding(setY.size());

						for (u64 itemIdx = 0; itemIdx < simple.mBins[binIdx].mBinRealSizes; ++itemIdx) //compute many F(k,xi)
							senderOprf.encode(binIdx*theirMaxBinSize + itemTheirIdx
								, &simple.mBins[binIdx].items[itemIdx], (u8*)&sendEncoding[itemIdx], sizeof(block));
					
						setY.emplace_back(mPrng.get<block>()); //add randome point
						sendEncoding.emplace_back(mPrng.get<block>());
						
						//poly
						std::cout << IoStream::lock;
						poly.getBlkCoefficients(simple.mMaxBinSize+1,sendEncoding, setY, coeffs);
						std::cout << IoStream::unlock;

						for (u64 c = 0; c < coeffs.size(); ++c)
							memcpy(sendBuff.data() + (k*itemTheirIdx*(simple.mMaxBinSize + 1)+c)* polyMaskBytes, (u8*)&coeffs[c], polyMaskBytes);
						
						//std::cout << IoStream::lock <<"r "<< binIdx << "\t" << itemTheirIdx << std::endl << IoStream::unlock;

					}

				}
				chl.asyncSend(std::move(sendBuff)); //done with sending P(x)

#if 0 //PEQT
				senderOprf.recvCorrection(chl, curStepSize*theirMaxBinSize);
				std::vector<block> sendEncoding(curStepSize*theirMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{
						senderOprf.encode(numOTs+binIdx*theirMaxBinSize + itemTheirIdx
							, &Ss[binIdx][itemTheirIdx], (u8*)&sendEncoding[k*simple.mMaxBinSize + itemTheirIdx], sizeof(block));
					}
				}

				/*u64 maskPEQTlength = mPsiSecParam / 8;
				std::vector<u8> recvBuff;
				chl.recv(recvBuff);
				if (recvBuff.size() != curStepSize*theirMaxBinSize*maskPEQTlength)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}*/

#endif
			}
		};



		//for (u64 i = 0; i <10; i += stepSize)
		//{

		//	auto curStepSize = std::min(stepSize, numOTs - i);
		//	senderOprf.recvCorrection(recvChl, curStepSize);

		//	std::vector<u8> sendBuff(curStepSize*(polyDegree)* poly.mNumBytes);

		//	for (u64 k = 0; k < curStepSize; ++k)
		//	{
		//		std::vector<block> setY(myInputSize + 1, Ss[i + k]);
		//		std::vector<block>sendEncoding(setY.size());

		//		for (u64 j = 0; j < myInputSize; ++j)
		//			senderOprf.encode(i + k, &inputs[j], (u8*)&sendEncoding[j], sizeof(block));

		//		//Poly
		//		polyNTL poly;
		//		poly.NtlPolyInit(polyNumBytes);
		//		setY[myInputSize] = mPrng.get<block>(); //add randome point
		//		sendEncoding[myInputSize] = mPrng.get<block>();

		//		poly.getBlkCoefficients(polyDegree, sendEncoding, setY, coeffs);

		//		for (u64 c = 0; c < coeffs.size(); ++c)
		//			memcpy(sendBuff.data() + (k*polyDegree + c)* poly.mNumBytes, (u8*)&coeffs[c], poly.mNumBytes);

		//		if (i + k == 0 || i + k == 1)
		//		{
		//			std::cout << "Ss[" << i + k << "]" << Ss[i + k] << "\t";
		//			std::cout << "coeffs[0]" << coeffs[0] << std::endl;

		//		}
		//	}
		//	recvChl.asyncSend(std::move(sendBuff));
		//}

		std::vector<std::thread> thrds(chls.size());
		for (u64 i = 0; i < thrds.size(); ++i)
		{
			thrds[i] = std::thread([=] {
				routine(i);
			});
		}

		for (auto& thrd : thrds)
			thrd.join();
	}
}
