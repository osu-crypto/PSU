#pragma once
// This file and the associated implementation has been placed in the public domain, waiving all copyright. No restrictions are placed on its use. 
#include <array>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Common/Timer.h>
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtSender.h"
#include "libPoly/polyNTL.h"

namespace osuCrypto
{

    class KrtwReceiver : public TimerAdapter
    {
    public:
     
		KrtwReceiver::KrtwReceiver()
		{
		}

		KrtwReceiver::~KrtwReceiver()
		{
		}

		bool mHasBase;

		u64 mNumOTs, mPolyNumBytes, mPolyDegree, mPsiSecParam;
		std::vector<block> mS;
		KkrtNcoOtSender senderOprf;
		polyNTL poly;
		PRNG mPrng;

		std::vector<block> mBaseOTRecv;
		BitVector mBaseChoice;

		void init(u64 psiSecParam, PRNG& prng, span<block> inputs, span<Channel> chls);
		void output(span<block> inputs, span<Channel> chls);

    };

}
