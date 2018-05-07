#include <iostream>

//using namespace std;
#include "tests_cryptoTools/UnitTests.h"
#include "libOTe_Tests/UnitTests.h"
#include <cryptoTools/gsl/span>

#include <cryptoTools/Common/Matrix.h>

#include <cryptoTools/Common/Defines.h>
using namespace osuCrypto;

#include "libOTe/TwoChooseOne/KosOtExtReceiver.h"
#include "libOTe/TwoChooseOne/KosOtExtSender.h"
#include "libOTe/TwoChooseOne/KosDotExtReceiver.h"
#include "libOTe/TwoChooseOne/KosDotExtSender.h"

#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Network/Endpoint.h>
#include <cryptoTools/Network/IOService.h>
#include <numeric>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Common/Log.h>


#include "libOTe/Tools/LinearCode.h"
#include "libOTe/Tools/bch511.h"
#include "libOTe/NChooseOne/Oos/OosNcoOtReceiver.h"
#include "libOTe/NChooseOne/Oos/OosNcoOtSender.h"
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtReceiver.h"
#include "libOTe/NChooseOne/Kkrt/KkrtNcoOtSender.h"

#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtSender.h"

#include "libOTe/NChooseK/AknOtReceiver.h"
#include "libOTe/NChooseK/AknOtSender.h"
#include "libOTe/TwoChooseOne/LzKosOtExtReceiver.h"
#include "libOTe/TwoChooseOne/LzKosOtExtSender.h"

#include "CLP.h"
#include "main.h"

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
#include "Poly/polyNTL.h"
#include "PsuDefines.h"

#include "PSU/KrtwSender.h"
#include "PSU/KrtwReceiver.h"
#include "Tools/SimpleIndex.h"

#include <thread>
#include <vector>

u64 disjontedSetSize = 10;
bool isTest=false;

template<typename ... Args>
std::string string_format(const std::string& format, Args ... args)
{
	size_t size = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
	std::unique_ptr<char[]> buf(new char[size]);
	std::snprintf(buf.get(), size, format.c_str(), args ...);
	return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

void Sender(u64 setSize, span<block> inputs, u64 numThreads)
{
	u64 psiSecParam = 40;

	PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
	PRNG prng1(_mm_set_epi32(4253465, 3434565, 234435, 23987025));

	// set up networking
	std::string name = "n";
	IOService ios;
	Endpoint ep1(ios, "localhost", 1212, EpMode::Server, name);

	std::vector<Channel> sendChls(numThreads);
	for (u64 i = 0; i < numThreads; ++i)
		sendChls[i] = ep1.addChannel("chl" + std::to_string(i), "chl" + std::to_string(i));

	KrtwSender sender;

	gTimer.reset();
	gTimer.setTimePoint("s start");
	sender.init(setSize, inputs.size(), 40, prng0, sendChls);
	gTimer.setTimePoint("s offline");

	/*std::cout << sender.mBaseOTSend[0][0] << "\t";
	std::cout << sender.mBaseOTSend[0][1] << "\n";
	std::cout << sender.mBaseOTRecv[0] << "\n";*/

	sender.output(inputs, sendChls);
	gTimer.setTimePoint("s finish");


	for (u64 i = 0; i < numThreads; ++i)
		sendChls[i].close();

	ep1.stop();	ios.stop();

}

void Receiver(u64 setSize, span<block> inputs,u64 numThreads)
{
	u64 psiSecParam = 40;

	PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
	PRNG prng1(_mm_set_epi32(4253465, 3434565, 234435, 23987025));

	// set up networking
	std::string name = "n";
	IOService ios;
	Endpoint ep0(ios, "localhost", 1212, EpMode::Client, name);

	std::vector<Channel> sendChls(numThreads), recvChls(numThreads);
	for (u64 i = 0; i < numThreads; ++i)
		recvChls[i] = ep0.addChannel("chl" + std::to_string(i), "chl" + std::to_string(i));

	KrtwReceiver recv;
	gTimer.reset();
	gTimer.setTimePoint("r start");
	recv.init(setSize, inputs.size(), 40, prng1,  recvChls);

		/*std::cout << recv.mBaseOTRecv[0] << "\n";
		std::cout << recv.mBaseOTSend[0][0] << "\t";
		std::cout << recv.mBaseOTSend[0][1] << "\n";*/

	recv.output(inputs, recvChls);
	gTimer.setTimePoint("r finish");
	std::cout << gTimer << std::endl;



	u64 dataSent = 0, dataRecv(0);
	for (u64 g = 0; g < recvChls.size(); ++g)
	{
		dataSent += recvChls[g].getTotalDataSent();
		dataRecv += recvChls[g].getTotalDataRecv();
		recvChls[g].resetStats();
	}

	std::cout << "Total Comm = " << string_format("%5.2f", (dataRecv + dataSent) / std::pow(2.0, 20)) << " MB\n";

	if (isTest)
	{
		std::cout << "recv.mDisjointedOutput.size(): " << recv.mDisjointedOutput.size() << std::endl;
		std::cout << "expectedDisjontedSetSize:      " << disjontedSetSize << std::endl;
	}

	for (u64 i = 0; i < numThreads; ++i)
		recvChls[i].close();

	ep0.stop(); ios.stop();

}

void PSU_Test_Impl()
{
	u64 setSize = 1 << 7, psiSecParam = 40, numThreads(1);

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
	std::cout << "intersection: " << sendSet[0] << "\n";
	std::cout << "intersection: " << sendSet[2] << "\n";


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
		recv.init(setSize, recvSet.size(), 40, prng1, recvChls);

		std::cout << recv.mBaseOTRecv[0] << "\n";

		std::cout << recv.mBaseOTSend[0][0] << "\t";
		std::cout << recv.mBaseOTSend[0][1] << "\n";

		recv.output(recvSet, recvChls);

	});

	sender.init(setSize, sendSet.size(), 40, prng0, sendChls);


	std::cout << sender.mBaseOTSend[0][0] << "\t";
	std::cout << sender.mBaseOTSend[0][1] << "\n";
	std::cout << sender.mBaseOTRecv[0] << "\n";

	sender.output(sendSet, sendChls);
	thrd.join();




	for (u64 i = 0; i < numThreads; ++i)
	{
		sendChls[i].close(); recvChls[i].close();
	}

	ep0.stop(); ep1.stop();	ios.stop();

}




void usage(const char* argv0)
{
	std::cout << "Error! Please use:" << std::endl;
	std::cout << "\t 1. For unit test: " << argv0 << " -t" << std::endl;
	std::cout << "\t 2. For simulation (2 terminal): " << std::endl;;
	std::cout << "\t\t Sender terminal: " << argv0 << " -r 0" << std::endl;
	std::cout << "\t\t Receiver terminal: " << argv0 << " -r 1" << std::endl;
}

int main(int argc, char** argv)
{
	
	u64 setSize = 1 << 12, numThreads = 1;


	if (argv[3][0] == '-' && argv[3][1] == 'n'
		&& argv[5][0] == '-' && argv[5][1] == 't')
	{
		setSize = 1 << atoi(argv[4]);
		numThreads = atoi(argv[6]);
	}


	std::cout << "SetSize: " << setSize << " vs " << setSize << "   |  numThreads: " << numThreads << "\t";

	PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
	PRNG prng1(_mm_set_epi32(4253233465, 334565, 0, 235));

	std::vector<block> sendSet(setSize), recvSet(setSize);
	for (u64 i = 0; i < setSize; ++i)
	{
		sendSet[i] = prng0.get<block>();
		recvSet[i] = sendSet[i];
	}

	for (u64 i = 0; i < disjontedSetSize; ++i)
		sendSet[i] = prng0.get<block>();

	//std::random_shuffle(sendSet.begin(), sendSet.begin(), prng0);


#if 0
	isTest = true;
	std::thread thrd = std::thread([&]() {
		Sender(setSize, sendSet, numThreads);
	});

	Receiver(setSize, recvSet, numThreads);

	thrd.join();
	return 0;
#endif

	if (argv[1][0] == '-' && argv[1][1] == 't') {
		
		isTest = true;
		std::thread thrd = std::thread([&]() {
			Sender(setSize,sendSet, numThreads);
		});

		Receiver(setSize, recvSet, numThreads);

		thrd.join();

	}
	else if (argv[1][0] == '-' && argv[1][1] == 'r' && atoi(argv[2]) == 0) {

		Sender(setSize, sendSet, numThreads);
	}
	else if (argv[1][0] == '-' && argv[1][1] == 'r' && atoi(argv[2]) == 1) {
		
		Receiver(setSize, recvSet, numThreads);
	}
	else {
		usage(argv[0]);
	}


	
  
	return 0;
}
