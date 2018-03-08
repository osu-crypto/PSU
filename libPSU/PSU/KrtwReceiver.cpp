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
	void KrtwReceiver::init(u64 psiSecParam, PRNG & prng, span<block> inputs, span<Channel> chls)
	{
		mPolyNumBytes = polyNumBytes;
		mPolyDegree = polyDegree;
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
		
		u64 numOTs = inputs.size(); //their input size
		u64 myInputSize = inputs.size();
		
		std::vector<block> coeffs, Ss(numOTs);
		for (u64 i = 0; i < numOTs; i++)
			Ss[i] = mPrng.get<block>();

		auto& recvChl = chls[0];

		senderOprf.init(numOTs, mPrng, recvChl);



		polyNTL poly;
		poly.NtlPolyInit(polyNumBytes);

		for (u64 i = 0; i <numOTs; i += stepSize)
		{

			auto curStepSize = std::min(stepSize, numOTs - i);
			senderOprf.recvCorrection(recvChl, curStepSize);

			std::vector<u8> sendBuff(curStepSize*(polyDegree)* poly.mNumBytes);

			for (u64 k = 0; k < curStepSize; ++k)
			{
				std::vector<block> setY(myInputSize + 1, Ss[i + k]);
				std::vector<block>sendEncoding(setY.size());

				for (u64 j = 0; j < myInputSize; ++j)
					senderOprf.encode(i + k, &inputs[j], (u8*)&sendEncoding[j], sizeof(block));

				//Poly
				setY[myInputSize] = mPrng.get<block>(); //add randome point
				sendEncoding[myInputSize] = mPrng.get<block>();

				poly.getBlkCoefficients(polyDegree, sendEncoding, setY, coeffs);

				for (u64 c = 0; c < coeffs.size(); ++c)
					memcpy(sendBuff.data() + (k*polyDegree + c)* poly.mNumBytes, (u8*)&coeffs[c], poly.mNumBytes);

				if (i + k == 0 || i + k == 1)
				{
					std::cout << "Ss[" << i + k << "]" << Ss[i + k] << "\t";
					std::cout << "coeffs[0]" << coeffs[0] << std::endl;

				}
			}
			recvChl.asyncSend(std::move(sendBuff));
		}
	}
}
