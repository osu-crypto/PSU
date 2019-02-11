#include "KrtwReceiver.h"

#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Crypto/Commit.h>
#include <cryptoTools/Network/Channel.h>
#include "libPSU/PsuDefines.h"

using namespace std;

namespace osuCrypto
{
	void KrtwReceiver::init(u64 myInputSize, u64 theirInputSize, u64 psiSecParam, PRNG & prng, span<Channel> chls)
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

		simple.init(myInputSize);
		Ss.resize(simple.mNumBins);

		theirMaxBinSize = simple.mMaxBinSize; //assume same set size, sender has mMaxBinSize, receiver has mMaxBinSize+1

		polyMaskBytes = (mPsiSecParam + log2(pow(simple.mMaxBinSize, 2)*simple.mNumBins) + 7) / 8;

		for (u64 i = 0; i < simple.mNumBins; i++)
		{
			Ss[i].resize(theirMaxBinSize);
			for (u64 j = 0; j < theirMaxBinSize; j++)
			{
				Ss[i][j] = ZeroBlock;
				block tem = mPrng.get<block>();
				memcpy((u8*)&Ss[i][j],(u8*)&tem, polyMaskBytes);
			}
		}

	/*	for (u64 j = 0; j < 6; j++)
			std::cout << "Ss [3][" << j << "]: " << Ss[3][j]  << "\n";
*/



	}
	void KrtwReceiver::output(span<block> inputs, span<Channel> chls)
	{
		u64 numThreads(chls.size());
		const bool isMultiThreaded = numThreads > 1;
		std::mutex mtx;

		//simple.print();

		//std::cout << "Receiver: " << simple.mMaxBinSize << "\t " <<simple.mNumBins<< std::endl ;

		u64	numOTs = simple.mNumBins*(theirMaxBinSize);

		std::vector<block> coeffs;

#ifdef DEBUG
		std::cout << IoStream::lock << "mBins[1].items[1]  " << simple.mBins[1].items[1] << std::endl << IoStream::unlock;
		std::cout << IoStream::lock << "Ss[1] " << Ss[1][1]<< std::endl << IoStream::unlock;
#endif
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

		gTimer.setTimePoint("r offline");

		simple.insertItems(inputs, numThreads);


#ifdef DEBUG
		std::cout << IoStream::lock << recvOTMsg[0] << std::endl << IoStream::unlock;
#endif
		//poly

		auto routine = [&](u64 t)
		{
			auto& chl = chls[t];
			u64 binStartIdx = simple.mNumBins * t / numThreads;
			u64 tempBinEndIdx = (simple.mNumBins * (t + 1) / numThreads);
			u64 binEndIdx = std::min(tempBinEndIdx, simple.mNumBins);
			
#ifdef _MSC_VER
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

				sendOprf.recvCorrection(chl, curStepSize*theirMaxBinSize); //OPRF

				std::vector<u8> sendBuff(curStepSize*theirMaxBinSize*(simple.mMaxBinSize)*polyMaskBytes);

				u64 iterSend = 0;

				//==========================PMT==========================
				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{

						//std::vector<block> setY(simple.mBins[binIdx].mBinRealSizes, Ss[binIdx][itemTheirIdx]);
						std::vector<block>sendEncoding(simple.mBins[binIdx].mBinRealSizes+1);
						u64 idxBot = simple.mBins[binIdx].mBinRealSizes;

						for (u64 itemIdx = 0; itemIdx < simple.mBins[binIdx].mBinRealSizes; itemIdx++) //compute many F(k,xi)
						{
							sendOprf.encode(binIdx*theirMaxBinSize + itemTheirIdx
								, &simple.mBins[binIdx].items[itemIdx], (u8*)&sendEncoding[itemIdx], sizeof(block));
						}


						//############ Global Item \bot ####################
						sendOprf.encode(binIdx*theirMaxBinSize + itemTheirIdx
							, &AllOneBlock, (u8*)&sendEncoding[idxBot], sizeof(block));


						//poly
#ifdef _MSC_VER
						std::cout << IoStream::lock;
						poly.getBlkCoefficients(simple.mMaxBinSize-1, sendEncoding, Ss[binIdx][itemTheirIdx], coeffs);
						std::cout << IoStream::unlock;
						/*coeffs.resize(simple.mMaxBinSize);
						for (u64 c = 0; c < coeffs.size(); ++c)
							coeffs[c] = ZeroBlock;*/

#else 
						poly.getBlkCoefficients(simple.mMaxBinSize-1, sendEncoding, Ss[binIdx][itemTheirIdx], coeffs);

#endif

						for (u64 c = 0; c < coeffs.size(); ++c)
						{
							memcpy(sendBuff.data() + iterSend, (u8*)&coeffs[c], polyMaskBytes);
							iterSend += polyMaskBytes;
						}

						//std::cout << IoStream::lock <<"r "<< binIdx << "\t" << itemTheirIdx << std::endl << IoStream::unlock;
					}
				}
				chl.asyncSend(std::move(sendBuff)); //done with sending P(x)
				//gTimer.setTimePoint("r sending coeffs");


				//==========================PEQT==========================
#if 0
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



				u64 maskPEQTlength = polyMaskBytes;// mPsiSecParam / 8;
				std::vector<u8> recvBuff;
				chl.recv(recvBuff); //receive OPRF(s*)
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
						
#ifdef DEBUG
						if (binIdx == 1 && itemTheirIdx == 1)
							std::cout << IoStream::lock << "recvS " << rcv << std::endl << IoStream::unlock;

#endif
						if (!memcmp(&rcv, &recvEncoding[k*theirMaxBinSize + itemTheirIdx], maskPEQTlength))  //get bit from PMT
						{
							bitPSU[k*theirMaxBinSize + itemTheirIdx] = 1;
							//std::cout << binIdx << "\t" << itemTheirIdx << "\t" << std::endl;
						}

						

					}
				}
				//gTimer.setTimePoint("r PMT output");


#endif


				//==========================receive S* directly==========================
#if 1
		
				u64 maskPEQTlength = polyMaskBytes;// mPsiSecParam / 8;
				std::vector<u8> recvBuff;
				chl.recv(recvBuff); //receive OPRF(s*)
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

#ifdef DEBUG
						if (binIdx == 1 && itemTheirIdx == 1)
							std::cout << IoStream::lock << "recvS " << rcv << std::endl << IoStream::unlock;

#endif
						if (!memcmp(&rcv, &Ss[binIdx][itemTheirIdx], maskPEQTlength))  //get bit from PMT
						{
							bitPSU[k*theirMaxBinSize + itemTheirIdx] = 1;
							//std::cout << binIdx << "\t" << itemTheirIdx << "\t" << std::endl;
						}

					}
				}
				//gTimer.setTimePoint("r PMT output");


#endif


				//==========================Rabin OT==========================
				sendBuff.resize(curStepSize*theirMaxBinSize*sizeof(u8));

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{
						u8 isOtMsgSwap = bitPSU[k*theirMaxBinSize + itemTheirIdx] ^ choicesOT[binIdx*theirMaxBinSize + itemTheirIdx];
							memcpy(sendBuff.data() + (k*theirMaxBinSize+ itemTheirIdx)
							, (u8*)&isOtMsgSwap, sizeof(u8));
					
#ifdef DEBUG
							if (binIdx == 1 && itemTheirIdx == 1)
								std::cout << IoStream::lock << bitPSU[k*theirMaxBinSize + itemTheirIdx] 
											<< "\t "<< choicesOT[binIdx*theirMaxBinSize + itemTheirIdx]
											<<"\t" <<unsigned(isOtMsgSwap)<< std::endl << IoStream::unlock;

#endif
					}
				}

				chl.asyncSend(std::move(sendBuff)); //done with sending choice OT

				//gTimer.setTimePoint("r Rabin OT");


#if 1 //PEQT

				//==========================compute PSU==========================

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

							psuItem = psuItem^recvOTMsg[binIdx*theirMaxBinSize + itemTheirIdx];

						/*	std::cout << "psuItem: " << psuItem << std::endl;
							std::cout << "itemTheirIdx: " << itemTheirIdx << std::endl;
							std::cout << "binIdx: " << binIdx << std::endl;
*/
							if (isMultiThreaded)
							{
								std::lock_guard<std::mutex> lock(mtx);
								{
									mDisjointedOutput.emplace_back(psuItem);
									/*std::cout << "itemTheirIdx: " << itemTheirIdx << std::endl;
									std::cout << "binIdx: " << binIdx << std::endl;*/
								}
							}
							else
							{
								mDisjointedOutput.emplace_back(psuItem);
							}
						}

					}
				}
#endif	

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

	void KrtwReceiver::outputNoOT(span<block> inputs, span<Channel> chls)
	{
		u64 numThreads(chls.size());
		const bool isMultiThreaded = numThreads > 1;
		std::mutex mtx;

		//simple.print();

		//std::cout << "Receiver: " << simple.mMaxBinSize << "\t " <<simple.mNumBins<< std::endl ;

		u64	numOTs = simple.mNumBins*(theirMaxBinSize);

		std::vector<std::vector<std::vector<block>>> coeffs(simple.mNumBins); //store p(.) for further brute force
		


#ifdef DEBUG
		std::cout << IoStream::lock << "mBins[1].items[1]  " << simple.mBins[1].items[1] << std::endl << IoStream::unlock;
		std::cout << IoStream::lock << "Ss[1] " << Ss[1][1] << std::endl << IoStream::unlock;
#endif
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

		gTimer.setTimePoint("r offline");

		simple.insertItems(inputs, numThreads);


#ifdef DEBUG
		std::cout << IoStream::lock << recvOTMsg[0] << std::endl << IoStream::unlock;
#endif
		//poly

		auto routine = [&](u64 t)
		{
			auto& chl = chls[t];
			u64 binStartIdx = simple.mNumBins * t / numThreads;
			u64 tempBinEndIdx = (simple.mNumBins * (t + 1) / numThreads);
			u64 binEndIdx = std::min(tempBinEndIdx, simple.mNumBins);

#ifdef _MSC_VER
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

				sendOprf.recvCorrection(chl, curStepSize*theirMaxBinSize); //OPRF

				std::vector<u8> sendBuff(curStepSize*theirMaxBinSize*(simple.mMaxBinSize)*polyMaskBytes);
				std::vector<std::vector<NTL::GF2EX>> p(curStepSize);

				u64 iterSend = 0;

				//==========================PMT==========================
				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					coeffs[binIdx].resize(theirMaxBinSize);
					p[k].resize(theirMaxBinSize);


					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{

						//std::vector<block> setY(simple.mBins[binIdx].mBinRealSizes, Ss[binIdx][itemTheirIdx]);
						std::vector<block>sendEncoding(simple.mBins[binIdx].mBinRealSizes + 1);
						u64 idxBot = simple.mBins[binIdx].mBinRealSizes;

						for (u64 itemIdx = 0; itemIdx < simple.mBins[binIdx].mBinRealSizes; itemIdx++) //compute many F(k,xi)
						{
							sendOprf.encode(binIdx*theirMaxBinSize + itemTheirIdx
								, &simple.mBins[binIdx].items[itemIdx], (u8*)&sendEncoding[itemIdx], sizeof(block));
						}


						//############ Global Item \bot ####################
						sendOprf.encode(binIdx*theirMaxBinSize + itemTheirIdx
							, &AllOneBlock, (u8*)&sendEncoding[idxBot], sizeof(block));


						//poly
#ifdef _MSC_VER
						std::cout << IoStream::lock;
						p[k][itemTheirIdx]=poly.getBlkCoefficients(simple.mMaxBinSize - 1, sendEncoding, Ss[binIdx][itemTheirIdx], coeffs[binIdx][itemTheirIdx]);
						std::cout << IoStream::unlock;
						/*coeffs.resize(simple.mMaxBinSize);
						for (u64 c = 0; c < coeffs.size(); ++c)
						coeffs[c] = ZeroBlock;*/

#else 
						p[k][itemTheirIdx] = poly.getBlkCoefficients(simple.mMaxBinSize - 1, sendEncoding, Ss[binIdx][itemTheirIdx], coeffs[binIdx][itemTheirIdx]);
#endif

						for (u64 c = 0; c < coeffs[binIdx][itemTheirIdx].size(); ++c)
						{
							memcpy(sendBuff.data() + iterSend, (u8*)&coeffs[binIdx][itemTheirIdx][c], polyMaskBytes);
							iterSend += polyMaskBytes;
						}

						//std::cout << IoStream::lock <<"r "<< binIdx << "\t" << itemTheirIdx << std::endl << IoStream::unlock;
					}
				}
				chl.asyncSend(std::move(sendBuff)); //done with sending P(x)
													//gTimer.setTimePoint("r sending coeffs");


				//==========================receive S* directly==========================
#if 1

				u64 maskSlength = polyMaskBytes;
				std::vector<u8> recvBuff;
				chl.recv(recvBuff); //receive OPRF(s*)
				if (recvBuff.size() != curStepSize*theirMaxBinSize*maskSlength)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemTheirIdx = 0; itemTheirIdx < theirMaxBinSize; ++itemTheirIdx)
					{
						block recvS; 
						memcpy((u8*)&recvS, recvBuff.data() + (k*theirMaxBinSize + itemTheirIdx)* maskSlength, maskSlength);
						
						if (memcmp(&recvS, &Ss[binIdx][itemTheirIdx], maskSlength))  //if s* \neq s => brute force x*
						{
				
							std::cout << binIdx << "\t" << itemTheirIdx << "\t"  << recvS << " - " << Ss[binIdx][itemTheirIdx] << std::endl;
							/*std::vector<block> blkRoots;
							poly.findRootsOffset(p[k][itemTheirIdx], recvS, blkRoots);*/

							block blkRoot;


							poly.findRootOffset(p[k][itemTheirIdx], recvS, blkRoot);

						}

					}
				}
				//gTimer.setTimePoint("r PMT output");


#endif




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
