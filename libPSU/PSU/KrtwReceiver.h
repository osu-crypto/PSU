#pragma once
// This file and the associated implementation has been placed in the public domain, waiving all copyright. No restrictions are placed on its use. 
#include <array>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Common/Timer.h>
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtSender.h"
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include "Poly/polyNTL.h"
#include "Tools/SimpleIndex.h"


namespace osuCrypto
{

    class KrtwReceiver : public TimerAdapter
    {
    public:
     
		
		bool mHasBase;

		u64 mNumOTs, mPolyNumBytes, mPolyDegree, mPsiSecParam;
		std::vector<block> mS;
		KkrtNcoOtSender sendOprf;
		KkrtNcoOtReceiver recvOprf;
		SimpleIndex simple;
		u64 polyMaskBytes;

		std::vector<std::vector<block>> Ss;
		u64 theirMaxBinSize;

		polyNTL poly;
		PRNG mPrng;

		std::vector<block> mBaseOTRecv;
		BitVector mBaseChoice;
		std::vector<std::array<block, 2>> mBaseOTSend;

		std::vector<block> mDisjointedOutput;

		void init(u64 myInputSize, u64 theirInputSize, u64 psiSecParam, PRNG& prng, span<Channel> chls);
		void output(span<block> inputs, span<Channel> chls);

    };

}
