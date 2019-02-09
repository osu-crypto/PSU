#include "KrtwSender.h"

#include <cryptoTools/Crypto/Commit.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/Timer.h>
#include "libOTe/Base/naor-pinkas.h"
#include "libPSU/PsuDefines.h"

namespace osuCrypto
{
    using namespace std;


	void KrtwSender::init(u64 myInputSize, u64 theirInputSize,u64 psiSecParam, PRNG & prng, span<Channel> chls)
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

		
		//std::cout << "baseCount "<< baseCount << std::endl;

		simple.init(myInputSize);
		theirMaxBinSize = simple.mMaxBinSize; //assume same set size, sender has mMaxBinSize, receiver has mMaxBinSize+1
		
		Sr.resize(simple.mNumBins);
		for (u64 i = 0; i < simple.mNumBins; i++)
			Sr[i].resize(simple.mMaxBinSize);

	}

	void KrtwSender::output(span<block> inputs, span<Channel> chls)
	{
		
		u64 numThreads(chls.size());
		//simple.print();
		//std::cout << IoStream::lock << "Sender: " << simple.mMaxBinSize << "\t " << simple.mNumBins<< std::endl << IoStream::unlock;

		u64	numOTs = simple.mNumBins*simple.mMaxBinSize;

		recvOprf.init( numOTs, mPrng, chls[0]); 
		sendOprf.init(numOTs, mPrng, chls[0]); //PEQT

		IknpOtExtSender sendIKNP;
		BitVector baseChoices(128);
		std::vector<block> baseRecv(128);

		baseChoices.copy(mBaseChoice, 0, 128);
		baseRecv.assign(mBaseOTRecv.begin(), mBaseOTRecv.begin() + 128);

		/*for (u64 i = 0; i < baseRecv.size(); i++)
		{
			baseChoices[i] = mBaseChoice[i];
			baseRecv[i] = mBaseOTRecv[i];
		}*/

		sendIKNP.setBaseOts(baseRecv, baseChoices);
		std::vector<std::array<block, 2>> sendOTMsg(numOTs);
		sendIKNP.send(sendOTMsg, mPrng, chls[0]);

		simple.insertItems(inputs, numThreads);
		//############ Global Item \bot ####################
		for (u64 i = 0; i < simple.mNumBins; ++i)
		{
			for (u64 j = simple.mBins[i].mBinRealSizes; j < simple.mMaxBinSize; ++j)
				simple.mBins[i].items.push_back(AllOneBlock);
		}
		



		//poly
		u64 polyMaskBytes = (mPsiSecParam + log2(pow(theirMaxBinSize,2)*simple.mNumBins) + 7) / 8;


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

				//==========================PMT==========================
				for (u64 k = 0; k < curStepSize; ++k) //OPRF
				{
					u64 binIdx = i + k;
					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						recvOprf.encode(binIdx*simple.mMaxBinSize + itemIdx
							, &simple.mBins[binIdx].items[itemIdx]
							, (u8*)&recvEncoding[k*simple.mMaxBinSize + itemIdx], sizeof(block));
					
						//if (binIdx == 3 && (itemIdx < 6))
						//{
						//	std::cout << itemIdx << ": " << simple.mBins[binIdx].items[itemIdx]
						//	<< "\t==rEnc==\t" << recvEncoding[k*simple.mMaxBinSize + itemIdx] << "\n";
						//}
					}
				}

				recvOprf.sendCorrection(chl, curStepSize*simple.mMaxBinSize);


#if 1 //poly

				std::vector<u8> recvBuff;
				chl.recv(recvBuff); //receive P(x)
				if (recvBuff.size() != curStepSize*simple.mMaxBinSize*(theirMaxBinSize)* polyMaskBytes)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}

#ifdef _MSC_VER
				std::cout << IoStream::lock;
				polyNTL poly;
				poly.NtlPolyInit(polyMaskBytes);//length=lambda +log(|Y|)
				std::cout << IoStream::unlock;
#else
				polyNTL poly;
				poly.NtlPolyInit(polyMaskBytes);//length=lambda +log(|Y|)
#endif

				u64 iterRecv = 0;

				std::vector<block> coeffs(theirMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{

						for (u64 c = 0; c < coeffs.size(); ++c)
						{
							memcpy((u8*)&coeffs[c], recvBuff.data() + iterRecv, polyMaskBytes);
							iterRecv += polyMaskBytes;
						}


#ifdef _MSC_VER
						std::cout << IoStream::lock;
						poly.evalPolynomial(coeffs, recvEncoding[k*simple.mMaxBinSize + itemIdx], Sr[binIdx][itemIdx]);
						std::cout << IoStream::unlock;
						//Sr[binIdx][itemIdx]=ZeroBlock;
#else
						poly.evalPolynomial(coeffs, recvEncoding[k*simple.mMaxBinSize + itemIdx], Sr[binIdx][itemIdx]);
#endif // _MSC_VER

					

#ifdef DEBUG
						if (binIdx == 3 && (itemIdx < 6))
						{
							std::cout << "Sr [3][" << itemIdx << "]: " << Sr[3][itemIdx] << "\t==rEnc\t"
								<< recvEncoding[k*simple.mMaxBinSize + itemIdx] << "\n";

						}
#endif
						
					}
			}
#endif
				//gTimer.setTimePoint("s compute s");


#if 0
				//==========================PEQT==========================

				sendOprf.recvCorrection(chl, curStepSize*simple.mMaxBinSize);
				std::vector<block> sendEncoding(curStepSize*simple.mMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						sendOprf.encode(binIdx*simple.mMaxBinSize + itemIdx
							, &Sr[binIdx][itemIdx], (u8*)&sendEncoding[k*simple.mMaxBinSize + itemIdx], sizeof(block));
#ifdef DEBUG
						if (binIdx == 1 && itemIdx == 1)
							std::cout << IoStream::lock << "sendEncoding " << sendEncoding[k*simple.mMaxBinSize + itemIdx] << std::endl << IoStream::unlock;
#endif
					}
				}


				u64 maskPEQTlength = mPsiSecParam / 8;
				std::vector<u8> sendBuff(curStepSize*simple.mMaxBinSize*maskPEQTlength);

				for (u64 c = 0; c < sendEncoding.size(); ++c)
					memcpy(sendBuff.data() + c* maskPEQTlength, (u8*)&sendEncoding[c], maskPEQTlength);

				chl.asyncSend(std::move(sendBuff)); //send OPRF(s*) == done with sending PEQT

				//gTimer.setTimePoint("s sending OPRF(s*)");
#endif

				//==========================send S* directly==========================
#if 1

				u64 maskPEQTlength = mPsiSecParam / 8;
				std::vector<u8> sendBuff(curStepSize*simple.mMaxBinSize*maskPEQTlength);

				int idxSendBuff = 0;

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						memcpy(sendBuff.data() + idxSendBuff, (u8*)&Sr[binIdx][itemIdx], maskPEQTlength);
						idxSendBuff += maskPEQTlength;

					}
				}

				chl.asyncSend(std::move(sendBuff)); //send (s*) 

#endif


#if 1 //Rabin OT


				//==========================Rabin OT==========================
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


						//
						if(itemIdx<simple.mBins[binIdx].mBinRealSizes)
							maskItem = simple.mBins[binIdx].items[itemIdx]^sendOTMsg[binIdx*simple.mMaxBinSize + itemIdx][isOtMsgSwap];
						else
							maskItem = AllOneBlock^sendOTMsg[binIdx*simple.mMaxBinSize + itemIdx][isOtMsgSwap];

						memcpy(sendBuff.data() + (k*simple.mMaxBinSize + itemIdx)* maskOTlength, (u8*)&maskItem, maskOTlength);

#ifdef DEBUG					
						if (binIdx == 1 && itemIdx == 1)
							std::cout << IoStream::lock << "isOtMsgSwap: " << unsigned(isOtMsgSwap) << std::endl << IoStream::unlock;

						if (binIdx == 1 && itemIdx == 1 ||
							binIdx == 3 && itemIdx == 23 ||
							binIdx == 4 && itemIdx == 19)
							std::cout << IoStream::lock << simple.mBins[binIdx].items[itemIdx]  <<" inter"<<std::endl << IoStream::unlock;
#endif
					}
				}
				chl.asyncSend(std::move(sendBuff)); //done with sending choice OT
				
				//gTimer.setTimePoint("s Rabin OT");

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


	void KrtwSender::outputNoOT(span<block> inputs, span<Channel> chls)
	{

		u64 numThreads(chls.size());
		//simple.print();
		//std::cout << IoStream::lock << "Sender: " << simple.mMaxBinSize << "\t " << simple.mNumBins<< std::endl << IoStream::unlock;

		u64	numOTs = simple.mNumBins*simple.mMaxBinSize;

		recvOprf.init(numOTs, mPrng, chls[0]);
		sendOprf.init(numOTs, mPrng, chls[0]); //PEQT

		IknpOtExtSender sendIKNP;
		BitVector baseChoices(128);
		std::vector<block> baseRecv(128);

		baseChoices.copy(mBaseChoice, 0, 128);
		baseRecv.assign(mBaseOTRecv.begin(), mBaseOTRecv.begin() + 128);

		/*for (u64 i = 0; i < baseRecv.size(); i++)
		{
		baseChoices[i] = mBaseChoice[i];
		baseRecv[i] = mBaseOTRecv[i];
		}*/

		sendIKNP.setBaseOts(baseRecv, baseChoices);
		std::vector<std::array<block, 2>> sendOTMsg(numOTs);
		sendIKNP.send(sendOTMsg, mPrng, chls[0]);

		simple.insertItems(inputs, numThreads);
		//############ Global Item \bot ####################
		for (u64 i = 0; i < simple.mNumBins; ++i)
		{
			for (u64 j = simple.mBins[i].mBinRealSizes; j < simple.mMaxBinSize; ++j)
				simple.mBins[i].items.push_back(AllOneBlock);
		}




		//poly
		u64 polyMaskBytes = (mPsiSecParam + log2(pow(theirMaxBinSize, 2)*simple.mNumBins) + 7) / 8;


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

				//==========================PMT==========================
				for (u64 k = 0; k < curStepSize; ++k) //OPRF
				{
					u64 binIdx = i + k;
					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						recvOprf.encode(binIdx*simple.mMaxBinSize + itemIdx
							, &simple.mBins[binIdx].items[itemIdx]
							, (u8*)&recvEncoding[k*simple.mMaxBinSize + itemIdx], sizeof(block));

						//if (binIdx == 3 && (itemIdx < 6))
						//{
						//	std::cout << itemIdx << ": " << simple.mBins[binIdx].items[itemIdx]
						//	<< "\t==rEnc==\t" << recvEncoding[k*simple.mMaxBinSize + itemIdx] << "\n";
						//}
					}
				}

				recvOprf.sendCorrection(chl, curStepSize*simple.mMaxBinSize);


#if 1 //poly

				std::vector<u8> recvBuff;
				chl.recv(recvBuff); //receive P(x)
				if (recvBuff.size() != curStepSize*simple.mMaxBinSize*(theirMaxBinSize)* polyMaskBytes)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}

#ifdef _MSC_VER
				std::cout << IoStream::lock;
				polyNTL poly;
				poly.NtlPolyInit(polyMaskBytes);//length=lambda +log(|Y|)
				std::cout << IoStream::unlock;
#else
				polyNTL poly;
				poly.NtlPolyInit(polyMaskBytes);//length=lambda +log(|Y|)
#endif

				u64 iterRecv = 0;

				std::vector<block> coeffs(theirMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{

						for (u64 c = 0; c < coeffs.size(); ++c)
						{
							memcpy((u8*)&coeffs[c], recvBuff.data() + iterRecv, polyMaskBytes);
							iterRecv += polyMaskBytes;
						}


#ifdef _MSC_VER
						std::cout << IoStream::lock;
						poly.evalPolynomial(coeffs, recvEncoding[k*simple.mMaxBinSize + itemIdx], Sr[binIdx][itemIdx]);
						std::cout << IoStream::unlock;
						//Sr[binIdx][itemIdx]=ZeroBlock;
#else
						poly.evalPolynomial(coeffs, recvEncoding[k*simple.mMaxBinSize + itemIdx], Sr[binIdx][itemIdx]);
#endif // _MSC_VER



#ifdef DEBUG
						if (binIdx == 3 && (itemIdx < 6))
						{
							std::cout << "Sr [3][" << itemIdx << "]: " << Sr[3][itemIdx] << "\t==rEnc\t"
								<< recvEncoding[k*simple.mMaxBinSize + itemIdx] << "\n";

						}
#endif

					}
				}
#endif
				//gTimer.setTimePoint("s compute s");


#if 0
				//==========================PEQT==========================

				sendOprf.recvCorrection(chl, curStepSize*simple.mMaxBinSize);
				std::vector<block> sendEncoding(curStepSize*simple.mMaxBinSize);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						sendOprf.encode(binIdx*simple.mMaxBinSize + itemIdx
							, &Sr[binIdx][itemIdx], (u8*)&sendEncoding[k*simple.mMaxBinSize + itemIdx], sizeof(block));
#ifdef DEBUG
						if (binIdx == 1 && itemIdx == 1)
							std::cout << IoStream::lock << "sendEncoding " << sendEncoding[k*simple.mMaxBinSize + itemIdx] << std::endl << IoStream::unlock;
#endif
					}
				}


				u64 maskPEQTlength = mPsiSecParam / 8;
				std::vector<u8> sendBuff(curStepSize*simple.mMaxBinSize*maskPEQTlength);

				for (u64 c = 0; c < sendEncoding.size(); ++c)
					memcpy(sendBuff.data() + c* maskPEQTlength, (u8*)&sendEncoding[c], maskPEQTlength);

				chl.asyncSend(std::move(sendBuff)); //send OPRF(s*) == done with sending PEQT

													//gTimer.setTimePoint("s sending OPRF(s*)");
#endif

			//==========================send S* directly==========================
#if 1

				u64 maskSlength = polyMaskBytes; //to brute force 
				std::vector<u8> sendBuff(curStepSize*simple.mMaxBinSize*maskSlength);

				int idxSendBuff = 0;

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;

					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
					{
						memcpy(sendBuff.data() + idxSendBuff, (u8*)&Sr[binIdx][itemIdx], maskSlength);
						idxSendBuff += maskSlength;

					}
				}

				chl.asyncSend(std::move(sendBuff)); //send (s*) 

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
