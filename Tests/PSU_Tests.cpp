#include "PSU_Tests.h"
#include "OT_Tests.h"

#include "libOTe/TwoChooseOne/OTExtInterface.h"

#include "libOTe/Tools/Tools.h"
#include "libOTe/Tools/LinearCode.h"
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Network/Endpoint.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Common/Log.h>

#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtSender.h"

#include "libOTe/TwoChooseOne/KosOtExtReceiver.h"
#include "libOTe/TwoChooseOne/KosOtExtSender.h"

#include "libOTe/TwoChooseOne/LzKosOtExtReceiver.h"
#include "libOTe/TwoChooseOne/LzKosOtExtSender.h"

#include "libOTe/TwoChooseOne/KosDotExtReceiver.h"
#include "libOTe/TwoChooseOne/KosDotExtSender.h"

#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtReceiver.h"
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtSender.h"
#include "libPoly/polyNTL.h"
#include "PsuDefines.h"

#include "PSU/KrtwSender.h"
#include "PSU/KrtwReceiver.h"
#include "Tools/SimpleIndex.h"

#include "Common.h"
#include <thread>
#include <vector>

#ifdef GetMessage
#undef GetMessage
#endif

#ifdef  _MSC_VER
#pragma warning(disable: 4800)
#endif //  _MSC_VER


using namespace osuCrypto;

namespace tests_libOTe
{

	inline void sse_trans(uint8_t *inp, int nrows, int ncols) {
#   define INP(x,y) inp[(x)*ncols/8 + (y)/8]
#   define OUT(x,y) inp[(y)*nrows/8 + (x)/8]
		int rr, cc, i, h;
		union { __m128i x; uint8_t b[16]; } tmp;
		__m128i vec;
		assert(nrows % 8 == 0 && ncols % 8 == 0);

		// Do the main body in 16x8 blocks:
		for (rr = 0; rr <= nrows - 16; rr += 16) {
			for (cc = 0; cc < ncols; cc += 8) {
				vec = _mm_set_epi8(
					INP(rr + 15, cc), INP(rr + 14, cc), INP(rr + 13, cc), INP(rr + 12, cc), INP(rr + 11, cc), INP(rr + 10, cc), INP(rr + 9, cc),
					INP(rr + 8, cc), INP(rr + 7, cc), INP(rr + 6, cc), INP(rr + 5, cc), INP(rr + 4, cc), INP(rr + 3, cc), INP(rr + 2, cc), INP(rr + 1, cc),
					INP(rr + 0, cc));
				for (i = 8; --i >= 0; vec = _mm_slli_epi64(vec, 1))
					*(uint16_t*)&OUT(rr, cc + i) = _mm_movemask_epi8(vec);
			}
		}
		if (rr == nrows) return;

		// The remainder is a block of 8x(16n+8) bits (n may be 0).
		//  Do a PAIR of 8x8 blocks in each step:
		for (cc = 0; cc <= ncols - 16; cc += 16) {
			vec = _mm_set_epi16(
				*(uint16_t const*)&INP(rr + 7, cc), *(uint16_t const*)&INP(rr + 6, cc),
				*(uint16_t const*)&INP(rr + 5, cc), *(uint16_t const*)&INP(rr + 4, cc),
				*(uint16_t const*)&INP(rr + 3, cc), *(uint16_t const*)&INP(rr + 2, cc),
				*(uint16_t const*)&INP(rr + 1, cc), *(uint16_t const*)&INP(rr + 0, cc));
			for (i = 8; --i >= 0; vec = _mm_slli_epi64(vec, 1)) {
				OUT(rr, cc + i) = h = _mm_movemask_epi8(vec);
				OUT(rr, cc + i + 8) = h >> 8;
			}
		}
		if (cc == ncols) return;

		//  Do the remaining 8x8 block:
		for (i = 0; i < 8; ++i)
			tmp.b[i] = INP(rr + i, cc);
		for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
			OUT(rr, cc + i) = _mm_movemask_epi8(tmp.x);
#undef INP
#undef OUT
	}
	


	void OT_Receive_Test(BitVector& choiceBits, gsl::span<block> recv, gsl::span<std::array<block, 2>>  sender)
	{

		for (u64 i = 0; i < choiceBits.size(); ++i)
		{

			u8 choice = choiceBits[i];
			const block & revcBlock = recv[i];
			//(i, choice, revcBlock);
			const block& senderBlock = sender[i][choice];

			//if (i%512==0) {
			//    std::cout << "[" << i << ",0]--" << sender[i][0] << std::endl;
			//    std::cout << "[" << i << ",1]--" << sender[i][1] << std::endl;
			//    std::cout << (int)choice << "-- " << recv[i] << std::endl;
			//}
			if (neq(revcBlock, senderBlock))
				throw UnitTestFail();

			if (eq(revcBlock, sender[i][1 ^ choice]))
				throw UnitTestFail();
		}

	}

    //\prod C(x)+C(xi)
	void Poly_IKNP_Test_Impl()
    {
        setThreadName("Sender");

        IOService ios(0);
        Endpoint ep0(ios, "127.0.0.1", 1212, EpMode::Server, "ep");
        Endpoint ep1(ios, "127.0.0.1", 1212, EpMode::Client, "ep");
        Channel senderChannel = ep1.addChannel("chl", "chl");
        Channel recvChannel = ep0.addChannel("chl", "chl");

        PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
        PRNG prng1(_mm_set_epi32(4253233465, 334565, 0, 235));

        u64 numOTs = 200;

        std::vector<block> recvMsg(numOTs), baseRecv(128);
        std::vector<std::array<block, 2>> sendMsg(numOTs), baseSend(128);
        BitVector choices(numOTs), baseChoice(128);
        choices.randomize(prng0);
        baseChoice.randomize(prng0);

        prng0.get((u8*)baseSend.data()->data(), sizeof(block) * 2 * baseSend.size());
        for (u64 i = 0; i < 128; ++i)
        {
            baseRecv[i] = baseSend[i][baseChoice[i]];
        }

        IknpOtExtSender sender;
        IknpOtExtReceiver recv;

        std::thread thrd = std::thread([&]() {
            recv.setBaseOts(baseSend);
            recv.receive(choices, recvMsg, prng0, recvChannel);
        });

        sender.setBaseOts(baseRecv, baseChoice);
        sender.send(sendMsg, prng1, senderChannel);
        thrd.join();

        //for (u64 i = 0; i < baseOTs.receiver_outputs.size(); ++i)
        //{
        //    std::cout << sender.GetMessage(i, 0) << " " << sender.GetMessage(i, 1) << "\n" << recv.GetMessage(1) << "  " << recv.mChoiceBits[i] << std::endl;
        //}
        OT_Receive_Test(choices, recvMsg, sendMsg);

		std::cout << "test\n";


		u64 beta = 32;

		static const u64 superBlkSize(8);
		u64 ell = superBlkSize*128;

		/*std::vector<std::array<block, superBlkSize>> codewords(128);

		for (u64 i = 0; i < 128; i++)
		{
		for (u64 j = 0; j < superBlkSize; j++)
		{
		prng1.get(&codewords[i][j], sizeof(block));
		}
		}*/
		
		block* x = new block[128];
		block* y = new block[128];

		prng1.get(y, 128);

		block a = AllOneBlock;
		a = a & (AllOneBlock^y[0]);
		a = a & (AllOneBlock^y[1]);
		a = a & (AllOneBlock^y[2]);
		a=a^ (y[0]&y[1]& y[2]);
		
		block b = AllOneBlock;
		b = b ^ y[0];
		b = b ^ y[1];
		b = b ^ y[2];
		b = b ^ (y[0] & y[1]);
		b = b ^ (y[0] & y[2]);
		b = b ^ (y[2] & y[1]);

		std::cout << "a= "<< a << std::endl;

		std::cout << b << std::endl;


		for (u64 i = 0; i < 128; i++)
		{
			x[i] = y[i];
		}
		x[1] = ZeroBlock;

		std::vector<std::array<block, superBlkSize>> xC(128), yC(128);

		std::array<block, superBlkSize> keys;
		prng1.get(keys.data(), keys.size());
		MultiKeyAES<superBlkSize> mMultiKeyAES;
		mMultiKeyAES.setKeys(keys);

		for (u64 i = 0; i < 128; i++)
		{
			std::array<block, superBlkSize> choice;
			for (u64 j = 0; j < superBlkSize; j++)
				choice[j] = x[i];
			mMultiKeyAES.ecbEncNBlocks(choice.data(), xC[i].data());

			for (u64 j = 0; j < superBlkSize; j++)
				prng1.get(&yC[i][j], 1);

				//choice[j] = y[i];
			//mMultiKeyAES.ecbEncNBlocks(choice.data(), yC[i].data());


		}
		
		std::cout << yC[0][0] << std::endl;
		std::cout << xC[0][0] << std::endl;

		std::cout << yC[0][1] << std::endl;
		std::cout << xC[0][1] << std::endl;

		std::cout << yC[1][1] << std::endl;
		std::cout << xC[1][1] << std::endl;

		std::vector<block> sBlock(superBlkSize);
		std::vector<block> prodC(superBlkSize);

		for (u64 i = 0; i < superBlkSize; i++)
		{
			sBlock[i] = AllOneBlock;
			prodC[i] = AllOneBlock;//C(y1)**C(y_n)
		}

		for (u64 j = 0; j < superBlkSize; j++)
		{
			for (u64 i = 0; i < 128; i++)
			{
				if(eq(yC[i][j],ZeroBlock))
					std::cout << i << "," <<j  << std::endl;
			}
		}

		for (u64 j = 0; j < superBlkSize; j++)
		{
			for (u64 i = 0; i < 128; i++)
			{
				sBlock[j] = sBlock[j] & (AllOneBlock^yC[i][j]);
				prodC[j] = prodC[j] &yC[i][j];
				std::cout << i << "," << j << ": " << yC[i][j] << std::endl;
				std::cout << i << "," << j << ": " <<  prodC[j] << std::endl;
			}
			sBlock[j] = sBlock[j] ^ prodC[j];
		}

		for (u64 j = 0; j < superBlkSize; j++)
		{
			std::cout << prodC[j] << std::endl;

		}

		
		std::cout << sBlock[0] << std::endl;

		BitVector s(ell);
		//s.fromBlock(&sBlock[0]);

		std::cout << s << std::endl;

		/*s.randomize(prng0);
		for (u64 i = 0; i < ell; i++)
			s[i] = 0;
		for (u64 i = 0; i < superBlkSize; i++)
			sBlock[i] = toBlock(s.data() + (i * sizeof(block)));*/

#if 1
		std::vector<block> q(ell);
		
		std::vector<std::array<block, 2>> t(ell);

		thrd = std::thread([&]() {
			recv.receive(s, q, prng0, recvChannel);
		});

		sender.send(t, prng1, senderChannel);
		thrd.join();

		OT_Receive_Test(s, q, t); //todo: get q & t directly from OT matrix
	
		

		std::vector<block> t0(ell);
		std::vector<block> t1(ell);

		for (u64 i = 0; i < ell; i++)
		{
			memcpy(&t0[i], &t[i][0], sizeof(block));
			memcpy(&t1[i], &t[i][1], sizeof(block));
		//	t0[i] = t[i][0];
		//	t1[i] = t[i][1];
		}
		std::cout << t0[0] << std::endl;
		std::cout << t[0][0] << std::endl;


		sse_trans((uint8_t*)&q[0], ell, 128);

		sse_trans((uint8_t*)&t0[0], ell, 128);
		sse_trans((uint8_t*)&t1[0], ell, 128);
	
		

		std::vector<std::array<block, superBlkSize>> corr(128);
		for (u64 i = 0; i < 128; i++)
		{
			for (u64 j = 0; j < superBlkSize; j++)
				corr[i][j] = xC[i][j] ^ t0[i*superBlkSize + j] ^ t1[i*superBlkSize + j];
		}

		for (u64 i = 0; i < 128; i++)
		{
			for (u64 j = 0; j < superBlkSize; j++)
				q[i*superBlkSize + j] = q[i*superBlkSize + j]^corr[i][j]& sBlock[j];
		}

		for (u64 i = 0; i < 128; i++)
		{
			for (u64 j = 0; j < superBlkSize; j++)
				q[i*superBlkSize + j] = q[i*superBlkSize + j] ^ prodC[j];
		}


		std::cout << t0[0] << std::endl;
		std::cout << q[0] << std::endl;

		std::cout << t0[1] << std::endl;
		std::cout << q[1] << std::endl;

		std::cout << t0[superBlkSize+1] << std::endl;
		std::cout << q[superBlkSize+1] << std::endl;


#endif

        senderChannel.close(); 
        recvChannel.close();
        ep1.stop();
        ep0.stop();
        ios.stop();
    }


	//PMT
	void PMT_Test_Impl1()
	{
		setThreadName("Sender"); 
		u64 sendSetSize = 20, recvSetSize = 20, psiSecParam = 40, maxBinSize = sendSetSize + 1;
		u64 polyNumBytes = 128 / 8, polyDegree = maxBinSize + 1;

		PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
		PRNG prng1(_mm_set_epi32(4253465, 3434565, 234435, 23987025));
		
		
		std::vector<block> sendSet(sendSetSize), recvSet(recvSetSize);
		for (u64 i = 0; i < sendSetSize; ++i)
			sendSet[i] = prng0.get<block>();

		for (u64 i = 0; i < recvSetSize; ++i)
			recvSet[i] = prng0.get<block>();

		sendSet[0] = recvSet[0];
		sendSet[2] = recvSet[2];
		

		// set up networking
		std::string name = "n";
		IOService ios;
		Session ep0(ios, "localhost", 1212, SessionMode::Server, name);
		Session ep1(ios, "localhost", 1212, SessionMode::Client, name);
		auto recvChl = ep1.addChannel(name, name);
		auto sendChl = ep0.addChannel(name, name);


#pragma region OPRF & poly
		// The total number that we wish to do
		u64 numOTs = recvSetSize;

		KkrtNcoOtSender senderOprf; senderOprf.configure(false, psiSecParam, 128);
		KkrtNcoOtReceiver recvOprf; recvOprf.configure(false, psiSecParam, 128);

		// the number of base OT that need to be done
		u64 baseCount = senderOprf.getBaseOTCount();

		// Fake some base OTs
		std::vector<block> baseRecv(baseCount);
		std::vector<std::array<block, 2>> baseSend(baseCount);
		BitVector baseChoice(baseCount);
		baseChoice.randomize(prng0);
		prng0.get((u8*)baseSend.data()->data(), sizeof(block) * 2 * baseSend.size());
		for (u64 i = 0; i < baseCount; ++i)
			baseRecv[i] = baseSend[i][baseChoice[i]];

		u64 stepSize = 5;
		std::vector<block> recvEncoding(numOTs), Ss(numOTs), Sr(numOTs);
		

		for (u64 i = 0; i < numOTs; i++)
			Ss[i]= prng0.get<block>();

		//recv
		auto thrd = std::thread([&]() {
			std::vector<block> coeffs;

			senderOprf.setBaseOts(baseRecv, baseChoice);
			senderOprf.init(numOTs, prng0, recvChl);
			polyNTL poly;
			poly.NtlPolyInit(polyNumBytes);

			u64 idx = 0;
			for (u64 i = 0; i < numOTs; i += stepSize)
			{
				auto curStepSize = std::min(stepSize, numOTs - i);
				senderOprf.recvCorrection(recvChl, curStepSize);


				std::vector<u8> sendBuff(curStepSize*(polyDegree)* poly.mNumBytes);

				for (u64 k = 0; k < curStepSize; ++k)
				{
					std::vector<block> setY(sendSetSize + 1, Ss[i + k]);
					std::vector<block>sendEncoding(setY.size());

					for (u64 j = 0; j < sendSetSize; ++j)
						senderOprf.encode(i + k, &sendSet[j], (u8*)&sendEncoding[j], sizeof(block));
					
					//Poly
					setY[sendSetSize] = prng0.get<block>(); //add randome point
					sendEncoding[sendSetSize] = prng0.get<block>(); 

					poly.getBlkCoefficients(polyDegree, sendEncoding, setY, coeffs);
					
					for (u64 c = 0; c < coeffs.size(); ++c)
						memcpy(sendBuff.data() + (k*polyDegree+c)* poly.mNumBytes, (u8*)&coeffs[c], poly.mNumBytes);

					if (i + k == 0 || i + k == 1)
					{
						std::cout << "Ss[" << i + k << "]" << Ss[i + k] << "\t";
						std::cout << "coeffs[0]" << coeffs[0] << std::endl;

					}
				}
				recvChl.asyncSend(std::move(sendBuff));
			}
		});
		
		recvOprf.setBaseOts(baseSend);
		recvOprf.init(numOTs, prng1, sendChl);
		polyNTL poly;
		poly.NtlPolyInit(polyNumBytes);

		for (u64 i = 0; i < numOTs; i += stepSize)
		{
			auto curStepSize = std::min(stepSize, numOTs - i);
			
			for (u64 k = 0; k < curStepSize; ++k)
				recvOprf.encode(i + k, &recvSet[i+k], (u8*)&recvEncoding[i+k], sizeof(block));
			
			recvOprf.sendCorrection(sendChl, curStepSize);

			//poly
			std::vector<u8> recvBuff;
			sendChl.recv(recvBuff);
			if (recvBuff.size() != curStepSize*(polyDegree)* poly.mNumBytes)
			{
				std::cout << "error @ " << (LOCATION) << std::endl;
				throw std::runtime_error(LOCATION);
			}

			std::vector<block> coeffs(polyDegree);

			for (u64 k = 0; k < curStepSize; ++k)
			{
				for (u64 c = 0; c < polyDegree; ++c)
					memcpy( (u8*)&coeffs[c], recvBuff.data() + (k*polyDegree + c)* poly.mNumBytes, poly.mNumBytes);
			
				poly.evalPolynomial(coeffs, recvEncoding[i + k], Sr[i + k]);
				if (i + k == 0 || i + k == 1)
				{
					std::cout << "Sr[" << i + k << "]" << Sr[i + k] << "\t";
					std::cout << "coeffs[0]" << coeffs[0] << std::endl;
				}
				

			}
		}

		thrd.join();
			
#pragma endregion




		sendChl.close(); recvChl.close();ep0.stop(); ep1.stop(); ios.stop();
	}


	void PMT_Test_Impl()
	{
		setThreadName("Sender");
		u64 setSize = 1 << 8, psiSecParam = 40, numThreads(1);

		PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
		PRNG prng1(_mm_set_epi32(4253465, 3434565, 234435, 23987025));


		std::vector<block> sendSet(setSize), recvSet(setSize);
		for (u64 i = 0; i < setSize; ++i)
		{
			sendSet[i] = prng0.get<block>();
			recvSet[i] = prng0.get<block>();
		}
		sendSet[0] = recvSet[0];
		sendSet[2] = recvSet[2];


		// set up networking
		std::string name = "n";
		IOService ios;
		Endpoint ep0(ios, "localhost", 1212, EpMode::Client, name);
		Endpoint ep1(ios, "localhost", 1212, EpMode::Server, name);

		std::vector<Channel> sendChls(numThreads), recvChls(numThreads);
		for (u64 i = 0; i < numThreads; ++i)
		{
			sendChls[i] = ep1.addChannel("chl" + std::to_string(i), "chl" + std::to_string(i));
			recvChls[i] = ep0.addChannel("chl" + std::to_string(i), "chl" + std::to_string(i));
		}


		KrtwSender sender;
		KrtwReceiver recv;
		auto thrd = std::thread([&]() {
			recv.init(40, prng1, recvSet, recvChls);

			std::cout << recv.mBaseOTRecv[0] << "\n";

			recv.output(recvSet, recvChls);

		});

		sender.init(40, prng0, sendSet, sendChls);
		

		std::cout << sender.mBaseOTSend[0][0] << "\n";
		std::cout << sender.mBaseOTSend[0][1] << "\n";

		sender.output(sendSet, sendChls);
		thrd.join();




		for (u64 i = 0; i < numThreads; ++i)
		{
			sendChls[i].close(); recvChls[i].close();
		}

		ep0.stop(); ep1.stop();	ios.stop();

	}

	void PSU_HashingParameters_Calculation() {
		SimpleIndex simpleIndex;
		std::vector<u64> logNumBalls{ 8, 12, 16, 20, 24 };
		std::vector<u64> lengthCodeWord{  424, 432, 440, 448, 448	};
		u64 statSecParam = 40, lengthItem=128, compSecParam=128;
		u64 commCost;
		double scale = 0,m=0;
		double iScaleStart = 0.03,iScaleEnd=0.08;

		for (u64 idxN = 0; idxN < logNumBalls.size(); idxN++)
		{
			u64 numBalls = 1 << logNumBalls[idxN];
			double iScale = iScaleStart;
			while (iScale < iScaleEnd)
			{
				u64 numBins = iScale*numBalls;
				u64 maxBinSize = simpleIndex.get_bin_size(numBins, numBalls, statSecParam);
				u64 curCommCost = numBins * (maxBinSize)*lengthCodeWord[idxN]
					+ numBins*maxBinSize * (maxBinSize + 2)*(statSecParam+log2(maxBinSize+2))
					+ numBins * maxBinSize*(lengthCodeWord[idxN] + statSecParam)
					+ numBins * maxBinSize*(compSecParam + lengthItem);
				

				if (iScale == iScaleStart)
				{
					commCost = curCommCost;
					scale = iScale;
					m = maxBinSize;
				}
				std::cout << iScale << "\t" << numBins << "\t" << maxBinSize <<"\t"
					<< curCommCost << "\t" << commCost << "\t";

				if (commCost > curCommCost)
				{
					commCost = curCommCost;
					scale = iScale;
					m = maxBinSize;

				}

				std::cout << scale << std::endl;

				//std::cout << iScale << "\t" << commCost <<"\t"<< curCommCost << std::endl;
				iScale += 0.0001;
			}
			std::cout << "##############" << std::endl;
			std::cout << logNumBalls[idxN] << "\t" << scale << "\t" << m << "\t"<< (commCost/8)*pow(10,-9)<< std::endl;
			std::cout << "##############" << std::endl;

		}
		
	}

	void Hashing_Test_Impl()
	{
		setThreadName("Sender");
		u64 setSize = 1<<8, psiSecParam = 40,  numThreads(2);

		PRNG prng(_mm_set_epi32(4253465, 3434565, 234435, 23987045));


		std::vector<block> set(setSize);
		for (u64 i = 0; i < set.size(); ++i)
			set[i] = prng.get<block>();

		SimpleIndex simple;
		simple.init(setSize);
		simple.insertItems(set,numThreads);
		simple.print();

	}

	void NTL_Poly_Test_Impl() {
		std::mutex mtx;

		auto routines = [&](u64 t)
		{
			
			polyNTL poly;
			poly.NtlPolyInit(128/8);
			PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));

			std::vector<block> setX(10);
			std::vector<block> setY(10, prng0.get<block>());

			block a= prng0.get<block>();
			for (u64 i = 0; i < 4; ++i)
			{
				setX[i] = prng0.get<block>();
			}

			block b = prng0.get<block>();
			for (u64 i = 5; i < setX.size(); ++i)
			{
				setX[i] = OneBlock;
			}

			setY[9] = prng0.get<block>();



			std::vector<block> coeffs;
			poly.getBlkCoefficients(11,setX, setY, coeffs);

			block y=ZeroBlock;
			poly.evalPolynomial(coeffs, setX[0],y);
			
			std::lock_guard<std::mutex> lock(mtx);
			std::cout << setY[0] << "\t" << y << std::endl;

		};

		std::vector<std::thread> thrds(1);
		for (u64 i = 0; i < thrds.size(); ++i)
		{
			thrds[i] = std::thread([=] {
				routines(i);
			});
		}

		for (auto& thrd : thrds)
			thrd.join();
	}

}