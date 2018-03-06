#include "KrtwSender.h"

#include <cryptoTools/Crypto/Commit.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/Timer.h>
#include "libOTe/Base/naor-pinkas.h"
#include "libPSU/PsuDefines.h"


namespace osuCrypto
{
    using namespace std;


	void KrtwSender::init(u64 psiSecParam, PRNG & prng, span<block> inputs, span<Channel> chls)
	{
		mPolyNumBytes = polyNumBytes;
		mPolyDegree = polyDegree;
		mPrng.SetSeed(prng.get<block>());
		recvOprf.configure(false, psiSecParam, 128);
		NaorPinkas baseOTs;
		mBaseOTSend.resize(recvOprf.getBaseOTCount());
		baseOTs.send(mBaseOTSend, mPrng, chls[0], 1);
		recvOprf.setBaseOts(mBaseOTSend);
	}

	void KrtwSender::output(span<block> inputs, span<Channel> chls)
	{
		
		u64 numOTs = inputs.size();
		
		auto& sendChl = chls[0];
		
		recvOprf.init(numOTs, mPrng, chls[0]);

		polyNTL poly;
		poly.NtlPolyInit(polyNumBytes);
		
		std::vector<block> recvEncoding(numOTs), Sr(numOTs);

		for (u64 i = 0; i < 10; i += stepSize)
		{
			auto curStepSize = std::min(stepSize, numOTs - i);

			for (u64 k = 0; k < curStepSize; ++k)
				recvOprf.encode(i + k, &inputs[i + k], (u8*)&recvEncoding[i + k], sizeof(block));

			recvOprf.sendCorrection(sendChl, curStepSize);

			//poly
#if 1
			std::vector<u8> recvBuff;
			sendChl.recv(recvBuff);
			if (recvBuff.size() != curStepSize*(mPolyDegree)* poly.mNumBytes)
			{
				std::cout << "error @ " << (LOCATION) << std::endl;
				throw std::runtime_error(LOCATION);
			}

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
	}

}
