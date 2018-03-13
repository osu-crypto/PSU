#include "KrtwSender.h"

#include <cryptoTools/Crypto/Commit.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/Timer.h>
#include "libOTe/Base/naor-pinkas.h"
#include "libPSU/PsuDefines.h"
#include "Tools/SimpleIndex.h"

namespace osuCrypto
{
    using namespace std;


	void KrtwSender::init(u64 psiSecParam, PRNG & prng, span<block> inputs, span<Channel> chls)
	{
		mPsiSecParam = psiSecParam;
		mPrng.SetSeed(prng.get<block>());
		recvOprf.configure(false, psiSecParam, 128);
		sendOprf.configure(false, psiSecParam, 128);

		u64 baseCount = sendOprf.getBaseOTCount();


		NaorPinkas baseOTs;
		mBaseOTSend.resize(baseCount);
		baseOTs.send(mBaseOTSend, mPrng, chls[0], 1);
		recvOprf.setBaseOts(mBaseOTSend);

		mBaseChoice.resize(baseCount);
		mBaseChoice.randomize(mPrng);
		mBaseOTRecv.resize(baseCount);
		baseOTs.receive(mBaseChoice, mBaseOTRecv, mPrng, chls[0], 1);
		sendOprf.setBaseOts(mBaseOTRecv, mBaseChoice);

		
		std::cout << "baseCount "<< baseCount << std::endl;

		

	}

	void KrtwSender::output(span<block> inputs, span<Channel> chls)
	{
		

		u64 numThreads(chls.size());
		
		SimpleIndex simple;
		simple.init(inputs.size(),false);
		simple.insertItems(inputs, numThreads);
		//simple.print();
		std::cout << IoStream::lock << "Sender: " << simple.mMaxBinSize << "\t " << simple.mNumBins
			<< std::endl << IoStream::unlock;




		u64 theirMaxBinSize = simple.mMaxBinSize + 1; //assume same set size, sender has mMaxBinSize, receiver has mMaxBinSize+1
		u64	numOTs = simple.mNumBins*simple.mMaxBinSize;
		
	
		std::vector<std::vector<block>> Sr(simple.mNumBins);
		for (u64 i = 0; i < simple.mNumBins; i++)
			Sr[i].resize(simple.mMaxBinSize);

		recvOprf.init( numOTs, mPrng, chls[0]); 
		sendOprf.init(numOTs, mPrng, chls[0]); //PEQT

		IknpOtExtSender sendIKNP;
		BitVector baseChoices(128);
		std::vector<block> baseRecv(128);

		for (u64 i = 0; i < baseRecv.size(); i++)
		{
			baseChoices[i] = mBaseChoice[i];
			baseRecv[i] = mBaseOTRecv[i];
		}

		sendIKNP.setBaseOts(baseRecv, baseChoices);
		std::vector<std::array<block, 2>> sendOTMsg(numOTs);
		sendIKNP.send(sendOTMsg, mPrng, chls[0]);

		std::cout << IoStream::lock << sendOTMsg[0][0] <<"\t" << sendOTMsg[0][1] << std::endl << IoStream::unlock;



		std::cout << IoStream::lock << "recvOprf.init done" << std::endl << IoStream::unlock;

		//poly
		u64 polyMaskBytes = 128 / 8;// (mPsiSecParam + log2(theirMaxBinSize + 1) + 7) / 8;


		auto routine = [&](u64 t)
		{
			auto& chl = chls[t];
			u64 binStartIdx = simple.mNumBins * t / numThreads;
			u64 tempBinEndIdx = (simple.mNumBins * (t + 1) / numThreads);
			u64 binEndIdx = std::min(tempBinEndIdx, simple.mNumBins);
			
			for (u64 i = binStartIdx; i < binEndIdx; i += stepSize)
			{
				auto curStepSize = std::min(stepSize, binEndIdx - i);
				std::vector<block> recvEncoding(curStepSize*simple.mMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemIdx = 0; itemIdx < simple.mBins[binIdx].mBinRealSizes; ++itemIdx)
					{
						recvOprf.encode(binIdx*simple.mMaxBinSize + itemIdx
							, &simple.mBins[binIdx].items[itemIdx]
							, (u8*)&recvEncoding[k*simple.mMaxBinSize + itemIdx], sizeof(block));
					}

					//TODO: send randome strings here
					for (u64 itemIdx = simple.mBins[binIdx].mBinRealSizes; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						block rnd = mPrng.get<block>();
						recvOprf.encode(binIdx*simple.mMaxBinSize + itemIdx
							, &simple.mBlkDefaut
							, (u8*)&recvEncoding[k*simple.mMaxBinSize + itemIdx], sizeof(block));
					}

				}

				recvOprf.sendCorrection(chl, curStepSize*simple.mMaxBinSize);


#if 1 //poly
				std::vector<u8> recvBuff;
				chl.recv(recvBuff);
				if (recvBuff.size() != curStepSize*simple.mMaxBinSize*(theirMaxBinSize + 1)* polyMaskBytes)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}

				std::cout << IoStream::lock;
				polyNTL poly;
				poly.NtlPolyInit(polyMaskBytes);//length=lambda +log(|Y|)
				std::cout << IoStream::unlock;

				std::vector<block> coeffs(theirMaxBinSize+1);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{

						for (u64 c = 0; c < coeffs.size(); ++c)
							memcpy((u8*)&coeffs[c], recvBuff.data() + (k*itemIdx*(theirMaxBinSize + 1) + c)* polyMaskBytes, polyMaskBytes);
							
						std::cout << IoStream::lock;
						poly.evalPolynomial(coeffs, recvEncoding[k*simple.mMaxBinSize + itemIdx], Sr[binIdx][itemIdx]);
						std::cout << IoStream::unlock;

						if (binIdx == 1 && itemIdx == 1)
						{
							std::cout << IoStream::lock << "mBins[1].items[1] " << simple.mBins[binIdx].items[itemIdx] << std::endl << IoStream::unlock;
							std::cout << IoStream::lock << "Sr[1] " << Sr[binIdx][itemIdx] << std::endl << IoStream::unlock;

						}
						//std::cout << IoStream::lock << "s " << binIdx << "\t"  << itemIdx << std::endl << IoStream::unlock;

					}

					//compute P(default) using same recvEnc
					/*for (u64 itemIdx = simple.mBins[binIdx].mBinRealSizes; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						Sr[binIdx][itemIdx] = Sr[binIdx][simple.mBins[binIdx].mBinRealSizes - 1];
					}*/


					/*if (binIdx == 1)
					{
						for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
						std::cout << IoStream::lock << "Sr[1].items[" << itemIdx << "] "
							<< "\t " << Sr[binIdx][itemIdx] << std::endl << IoStream::unlock;

					}*/

			}
#endif

#if 1 //PEQT


				sendOprf.recvCorrection(chl, curStepSize*simple.mMaxBinSize);
				std::vector<block> sendEncoding(curStepSize*simple.mMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						sendOprf.encode(binIdx*simple.mMaxBinSize + itemIdx
							, &Sr[binIdx][itemIdx], (u8*)&sendEncoding[k*simple.mMaxBinSize + itemIdx], sizeof(block));
					
						if (binIdx == 1 && itemIdx == 1)
						{
							std::cout << IoStream::lock << "sendEncoding " << sendEncoding[k*simple.mMaxBinSize + itemIdx] << std::endl << IoStream::unlock;

						}
					}
				}


				u64 maskPEQTlength = mPsiSecParam / 8;
				std::vector<u8> sendBuff(curStepSize*simple.mMaxBinSize*maskPEQTlength);

				for (u64 c = 0; c < sendEncoding.size(); ++c)
					memcpy(sendBuff.data() + c* maskPEQTlength, (u8*)&sendEncoding[c], maskPEQTlength);

				chl.asyncSend(std::move(sendBuff)); //done with sending PEQT


				

#endif
				chl.recv(recvBuff); //receive isOtMsgSwap
				if (recvBuff.size() != curStepSize*simple.mMaxBinSize)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}


				u64 maskOTlength = 128 / 8; //get final output
				sendBuff.resize(curStepSize*simple.mMaxBinSize*maskOTlength);
				block maskItem;

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						u8 isOtMsgSwap;
						memcpy((u8*)&isOtMsgSwap, recvBuff.data() + (k*simple.mMaxBinSize + itemIdx), sizeof(u8));
					
						if (binIdx == 1 && itemIdx == 1)
							std::cout << IoStream::lock << "isOtMsgSwap: "  << unsigned(isOtMsgSwap) << std::endl << IoStream::unlock;

						//
						if(itemIdx<simple.mBins[binIdx].mBinRealSizes)
							maskItem = simple.mBins[binIdx].items[itemIdx] + sendOTMsg[binIdx*simple.mMaxBinSize + itemIdx][isOtMsgSwap];
						else
							maskItem = simple.mBlkDefaut + sendOTMsg[binIdx*simple.mMaxBinSize + itemIdx][isOtMsgSwap];

						memcpy(sendBuff.data() + (k*simple.mMaxBinSize + itemIdx)* maskOTlength, (u8*)&maskItem, maskOTlength);

					
						if (binIdx == 1 && itemIdx == 1 ||
							binIdx == 3 && itemIdx == 23 ||
							binIdx == 4 && itemIdx == 19)
							std::cout << IoStream::lock << simple.mBins[binIdx].items[itemIdx]  <<" inter"<<std::endl << IoStream::unlock;

					}
				}
				chl.asyncSend(std::move(sendBuff)); //done with sending choice OT
			}


		};

	
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
