#pragma once
// This file and the associated implementation has been placed in the public domain, waiving all copyright. No restrictions are placed on its use.  
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Network/Channel.h>
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtReceiver.h"
#include "libPoly/polyNTL.h"

#include <array>
namespace osuCrypto {

	class KrtwSender :public TimerAdapter
	{
	public:

		KrtwSender::KrtwSender()
		{
		}

		KrtwSender::~KrtwSender()
		{
		}

		bool mHasBase;

		u64 mNumOTs, mPolyNumBytes, mPolyDegree, mStepSize, mPsiSecParam;
		std::vector<block> mS;
		KkrtNcoOtReceiver recvOprf; 
		polyNTL poly;
		PRNG mPrng;
		
		std::vector<std::array<block, 2>> mBaseOTSend;

		void init(u64 psiSecParam, PRNG& prng, span<block> inputs, span<Channel> chls);
		void output(span<block> inputs, span<Channel> chls);
	};
}

