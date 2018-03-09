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
		NaorPinkas baseOTs;
		mBaseOTSend.resize(recvOprf.getBaseOTCount());
		baseOTs.send(mBaseOTSend, mPrng, chls[0], 1);
		recvOprf.setBaseOts(mBaseOTSend);
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
		std::vector<block> recvEncoding(numOTs), Sr(numOTs);

		recvOprf.init(numOTs, mPrng, chls[0]);

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

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 binIdx = i + k;
					for (u64 itemIdx = 0; itemIdx < simple.mMaxBinSize; ++itemIdx)
						recvOprf.encode(binIdx*simple.mMaxBinSize+itemIdx
										, &simple.mBins[binIdx].items[itemIdx]
										, (u8*)&recvEncoding[binIdx*simple.mMaxBinSize + itemIdx], sizeof(block));

				}

				recvOprf.sendCorrection(chl, curStepSize*simple.mMaxBinSize);

				//poly

				std::vector<u8> recvBuff;
				chl.recv(recvBuff);
				if (recvBuff.size() != curStepSize*simple.mMaxBinSize*(theirMaxBinSize+1)* polyMaskBytes)
				{
					std::cout << "error @ " << (LOCATION) << std::endl;
					throw std::runtime_error(LOCATION);
				}
#if 0
				polyNTL poly;
				poly.NtlPolyInit(((mPsiSecParam + log2(simple.mMaxBinSize + 2) + 7) / 8));
				std::vector<block> coeffs(mPolyDegree);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					for (u64 c = 0; c < mPolyDegree; ++c)
						memcpy((u8*)&coeffs[c], recvBuff.data() + (k*mPolyDegree + c)* poly.mNumBytes, poly.mNumBytes);

					poly.evalPolynomial(coeffs, recvEncoding[i + k], Sr[i + k]);
					if (i + k == 0 || i + k == 1)
					{
						std::cout << "Sr[" << i + k << "]" << Sr[i + k] << "\t";
						std::cout << "coeffs[0]" << coeffs[0] << std::endl;
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

}
