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

		sendOprf.configure(false, psiSecParam, 128);
		recvOprf.configure(false, psiSecParam, 128);

		u64 baseCount= sendOprf.getBaseOTCount();

		mBaseChoice.resize(baseCount);
		mBaseChoice.randomize(mPrng);
		mBaseOTRecv.resize(baseCount);
		NaorPinkas baseOTs;
		baseOTs.receive(mBaseChoice, mBaseOTRecv, mPrng, chls[0], 1);
		sendOprf.setBaseOts(mBaseOTRecv, mBaseChoice);


		mBaseOTSend.resize(baseCount);
		baseOTs.send(mBaseOTSend, mPrng, chls[0], 1);
		recvOprf.setBaseOts(mBaseOTSend);

		

	}
	void KrtwReceiver::output(span<block> inputs, span<Channel> chls)
	{
		u64 numThreads(chls.size());
		const bool isMultiThreaded = numThreads > 1;

		std::mutex mtx;

		SimpleIndex simple;
		simple.init(inputs.size(),true);
		simple.insertItems(inputs, numThreads);
		//simple.print();

		//std::cout << "Receiver: " << simple.mMaxBinSize << "\t " <<simple.mNumBins<< std::endl ;

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

		sendOprf.init(numOTs, mPrng, chls[0]);
		recvOprf.init(numOTs, mPrng, chls[0]);//PEQT
	
		
		IknpOtExtReceiver recvIKNP;
		std::vector<std::array<block, 2>> baseOTSend(128);

		for (u64 i = 0; i < baseOTSend.size(); i++)
		{
			baseOTSend[i][0] = mBaseOTSend[i][0];
			baseOTSend[i][1] = mBaseOTSend[i][1];
		}
		recvIKNP.setBaseOts(baseOTSend);

		BitVector choicesOT(numOTs); choicesOT.randomize(mPrng);
		std::vector<block> recvOTMsg(numOTs);
		recvIKNP.receive(choicesOT, recvOTMsg, mPrng, chls[0]);

		std::cout << IoStream::lock << recvOTMsg[0] << std::endl << IoStream::unlock;

		//poly
		u64 polyMaskBytes = (mPsiSecParam + log2(simple.mMaxBinSize + 1) + 7) / 8;

		auto routine = [&](u64 t)
		{
			auto& chl = chls[t];
			u64 binStartIdx = simple.mNumBins * t / numThreads;
			u64 tempBinEndIdx = (simple.mNumBins * (t + 1) / numThreads);
			u64 binEndIdx = std::min(tempBinEndIdx, simple.mNumBins);
			
#ifdef NTL_Threads
			std::cout << IoStream::lock;
			polyNTL poly;
			poly.NtlPolyInit(polyMaskBytes);//length=lambda +log(|Y|)
			std::cout << IoStream::unlock;
#else
			polyNTL poly;
			poly.NtlPolyInit(polyMaskBytes);//length=lambda +log(|Y|)
#endif

			for (u64 i = binStartIdx; i < binEndIdx; i += stepSize)
			{
				auto curStepSize = std::min(stepSize, binEndIdx - i);

				sendOprf.recvCorrection(chl, curStepSize*theirMaxBinSize);

				std::vector<u8> sendBuff(curStepSize*theirMaxBinSize*(simple.mMaxBinSize + 1)*polyMaskBytes);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{

						std::vector<block> setY(simple.mBins[binIdx].mBinRealSizes, Ss[binIdx][itemTheirIdx]);
						std::vector<block>sendEncoding(setY.size());

						for (u64 itemIdx = 0; itemIdx < simple.mBins[binIdx].mBinRealSizes; ++itemIdx) //compute many F(k,xi)
							sendOprf.encode(binIdx*theirMaxBinSize + itemTheirIdx
								, &simple.mBins[binIdx].items[itemIdx], (u8*)&sendEncoding[itemIdx], sizeof(block));

						setY.emplace_back(mPrng.get<block>()); //add randome point
						sendEncoding.emplace_back(mPrng.get<block>());

						//poly
#ifdef NTL_Threads
						std::cout << IoStream::lock;
						poly.getBlkCoefficients(simple.mMaxBinSize + 1, sendEncoding, setY, coeffs);
						std::cout << IoStream::unlock;
#else
						poly.getBlkCoefficients(simple.mMaxBinSize + 1, sendEncoding, setY, coeffs);
#endif
						for (u64 c = 0; c < coeffs.size(); ++c)
							memcpy(sendBuff.data() + (k*itemTheirIdx*(simple.mMaxBinSize + 1) + c)* polyMaskBytes, (u8*)&coeffs[c], polyMaskBytes);

						//std::cout << IoStream::lock <<"r "<< binIdx << "\t" << itemTheirIdx << std::endl << IoStream::unlock;
					}
				}
				chl.asyncSend(std::move(sendBuff)); //done with sending P(x)

#if 1 //PEQT

				std::vector<block> recvEncoding(curStepSize*theirMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{
						//using second part of recvOprf (start at numOTs idx)
						recvOprf.encode(binIdx*theirMaxBinSize + itemTheirIdx
							, &Ss[binIdx][itemTheirIdx]
							, (u8*)&recvEncoding[k*theirMaxBinSize + itemTheirIdx], sizeof(block));
					
#ifdef DEBUG
						if (binIdx == 1 && itemTheirIdx == 1)
							std::cout << IoStream::lock << "recvEncoding " << recvEncoding[k*theirMaxBinSize + itemTheirIdx] << std::endl << IoStream::unlock;
#endif
					}
				}

				recvOprf.sendCorrection(chl, curStepSize*theirMaxBinSize);

				u64 maskPEQTlength = mPsiSecParam / 8;
				std::vector<u8> recvBuff;
				chl.recv(recvBuff);
				if (recvBuff.size() != curStepSize*theirMaxBinSize*maskPEQTlength)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}

				BitVector bitPSU(curStepSize*theirMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{
						block rcv;
						memcpy((u8*)&rcv, recvBuff.data() + (k*theirMaxBinSize + itemTheirIdx)* maskPEQTlength, maskPEQTlength);
						
						if (binIdx == 1 && itemTheirIdx == 1)
							std::cout << IoStream::lock << "rcv " << rcv << std::endl << IoStream::unlock;

						if (!memcmp(&rcv, &recvEncoding[k*theirMaxBinSize + itemTheirIdx], maskPEQTlength))
						{
							bitPSU[k*theirMaxBinSize + itemTheirIdx] = 1;
							//std::cout << binIdx << "\t" << itemTheirIdx << "\t" << std::endl;
						}

					}
				}

#endif

				sendBuff.resize(curStepSize*theirMaxBinSize*sizeof(u8));

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{
						u8 isOtMsgSwap = bitPSU[k*theirMaxBinSize + itemTheirIdx] ^ choicesOT[binIdx*theirMaxBinSize + itemTheirIdx];
							memcpy(sendBuff.data() + (k*theirMaxBinSize+ itemTheirIdx)
							, (u8*)&isOtMsgSwap, sizeof(u8));
					
							if (binIdx == 1 && itemTheirIdx == 1)
								std::cout << IoStream::lock << bitPSU[k*theirMaxBinSize + itemTheirIdx] 
											<< "\t "<< choicesOT[binIdx*theirMaxBinSize + itemTheirIdx]
											<<"\t" <<unsigned(isOtMsgSwap)<< std::endl << IoStream::unlock;

					
					}
				}



				chl.asyncSend(std::move(sendBuff)); //done with sending choice OT


				u64 maskOTlength = 128 / 8;
				chl.recv(recvBuff);
				if (recvBuff.size() != curStepSize*theirMaxBinSize*maskOTlength)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}

				
				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{
						if (!bitPSU[k*theirMaxBinSize + itemTheirIdx])
						{
							block psuItem;
							memcpy((u8*)&psuItem, recvBuff.data() + (k*theirMaxBinSize + itemTheirIdx)* maskOTlength, maskOTlength);

							psuItem = psuItem + recvOTMsg[binIdx*theirMaxBinSize + itemTheirIdx];


							if (isMultiThreaded)
							{
								std::lock_guard<std::mutex> lock(mtx);
								PsuOutput.emplace_back(psuItem);
							}
							else
							{
								PsuOutput.emplace_back(psuItem);
							}
						}

					}
				}
			

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

		std::cout << "PsuOutput.size() " <<PsuOutput.size() << std::endl;

	}
}
