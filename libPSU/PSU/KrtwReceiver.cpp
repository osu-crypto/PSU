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

		std::vector<block> coeffs, Ss(numOTs);
		for (u64 i = 0; i < numOTs; i++)
			Ss[i] = mPrng.get<block>();

		senderOprf.init(numOTs, mPrng, chls[0]);

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

			for (u64 i = binStartIdx; i < binEndIdx; i += stepSize)
			{
				auto curStepSize = std::min(stepSize, binEndIdx - i);
				std::vector<u8> sendBuff(curStepSize*theirMaxBinSize*(simple.mMaxBinSize + 1)*polyMaskBytes);

				senderOprf.recvCorrection(chl, curStepSize*theirMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{
						std::vector<block> setY(simple.mMaxBinSize+1, Ss[binIdx*theirMaxBinSize+ itemTheirIdx]);
						std::vector<block>sendEncoding(setY.size());

						for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx) //compute many F(k,xi)
							senderOprf.encode(binIdx*theirMaxBinSize + itemTheirIdx
								, &simple.mBins[binIdx].items[itemIdx], (u8*)&sendEncoding[itemIdx], sizeof(block));
					
					

						//poly
						for (size_t ii = 0; ii < setY.size(); ii++)
						{
							std::cout << IoStream::lock << sendEncoding[ii] << "\t" << setY[ii] << std::endl << IoStream::unlock;

							setY[ii] = mPrng.get<block>(); //add randome point
							sendEncoding[ii] = mPrng.get<block>();
						}
						setY[setY.size()-1] = mPrng.get<block>(); //add randome point
						sendEncoding[setY.size() - 1] = mPrng.get<block>();

						polyNTL poly;
						poly.NtlPolyInit(polyMaskBytes); //length=lambda +log(|Y|)
					
						poly.getBlkCoefficients(sendEncoding, setY, coeffs);

						for (u64 c = 0; c < coeffs.size(); ++c)
							memcpy(sendBuff.data() + (k*itemTheirIdx + c)* polyMaskBytes
								, (u8*)&coeffs[c], polyMaskBytes);
					}
				}
				
				chl.asyncSend(std::move(sendBuff));
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
