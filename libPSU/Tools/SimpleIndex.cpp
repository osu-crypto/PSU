#include "SimpleIndex.h"
#include "cryptoTools/Crypto/sha1.h"
#include "cryptoTools/Crypto/PRNG.h"
#include <random>
#include "cryptoTools/Common/Log.h"
#include "cryptoTools/Common/CuckooIndex.h"
#include <numeric>
//#include <boost/math/special_functions/binomial.hpp>
//#include <boost/multiprecision/cpp_bin_float.hpp>

namespace osuCrypto
{


    void SimpleIndex::print()
    {
		std::cout << "mMaxBinSize=" << mMaxBinSize << std::endl;
		std::cout << "mNumBins=" << mNumBins << std::endl;
        for (u64 i = 0; i < mBins.size(); ++i)
        {
            std::cout << "Bin #" << i << std::endl;

            std::cout << " contains " << mBins[i].mBinRealSizes << " elements" << std::endl;

            for (u64 j = 0; j < mBins[i].items.size(); ++j)
                std::cout << "    idx= " << j << ": " << mBins[i].items[j] << std::endl;
			
            std::cout << std::endl;
        }

        std::cout << std::endl;
    }

	void SimpleIndex::print(u64 numBin)
	{
		std::cout << "mMaxBinSize=" << mMaxBinSize << std::endl;
		std::cout << "mNumBins=" << mNumBins << std::endl;
		u64 i = numBin;
		{
			std::cout << "Bin #" << i << std::endl;

			std::cout << " contains " << mBins[i].mBinRealSizes << " elements" << std::endl;

			for (u64 j = 0; j < mBins[i].items.size(); ++j)
				std::cout << "    idx= " << j << ": " << mBins[i].items[j] << std::endl;


			std::cout << std::endl;
		}

		std::cout << std::endl;
	}

#if 0
   // template<unsigned int N = 16>
    double getBinOverflowProb(u64 numBins, u64 numBalls, u64 binSize, double epsilon = 0.0001)
    {
        if (numBalls <= binSize)
            return std::numeric_limits<double>::max();

        if (numBalls > std::numeric_limits<i32>::max())
        {
            auto msg = ("boost::math::binomial_coefficient(...) only supports " + std::to_string(sizeof(unsigned) * 8) + " bit inputs which was exceeded." LOCATION);
            std::cout << msg << std::endl;
            throw std::runtime_error(msg);
        }

        //std::cout << numBalls << " " << numBins << " " << binSize << std::endl;
        typedef boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<16>> T;
        T sum = 0.0;
        T sec = 0.0;// minSec + 1;
        T diff = 1;
        u64 i = binSize + 1;


        while (diff > T(epsilon) && numBalls >= i /*&& sec > minSec*/)
        {
            sum += numBins * boost::math::binomial_coefficient<T>(i32(numBalls), i32(i))
                * boost::multiprecision::pow(T(1.0) / numBins, i) * boost::multiprecision::pow(1 - T(1.0) / numBins, numBalls - i);

            //std::cout << "sum[" << i << "] " << sum << std::endl;

            T sec2 = boost::multiprecision::log2(sum);
            diff = boost::multiprecision::abs(sec - sec2);
            //std::cout << diff << std::endl;
            sec = sec2;

            i++;
        }

        return std::max<double>(0, (double)-sec);
    }


    u64 get_bin_size(u64 numBins, u64 numBalls, u64 statSecParam)
    {

        auto B = std::max<u64>(1, numBalls / numBins);

        double currentProb = getBinOverflowProb(numBins, numBalls, B);
        u64 step = 1;

        bool doubling = true;

        while (currentProb < statSecParam || step > 1)
        {
            if (!step)
                throw std::runtime_error(LOCATION);


            if (statSecParam > currentProb)
            {
                if (doubling) step = std::max<u64>(1, step * 2);
                else          step = std::max<u64>(1, step / 2);

                B += step;
            }
            else
            {
                doubling = false;
                step = std::max<u64>(1, step / 2);
                B -= step;
            }
            currentProb = getBinOverflowProb(numBins, numBalls, B);
        }

        return B;
    }

#endif

	u64 SimpleIndex::getMaxBinSize(u64 numBalls) {
		if (numBalls <= 1 << 8)
			return 63;
		else if (numBalls <= 1 << 10)
			return 58; //10	0.055	58	1.81555 Mb
		else if (numBalls <= 1 << 12)
			return 63; //12	0.05	63	7.84829 Mb
		else if (numBalls <= 1 << 14)
			return 62; //14	0.053	62	33.4393 Mb
		else if (numBalls <= 1 << 16)
			return 60; //16	0.058	60	141.783 Mb
		else if (numBalls <= 1 << 18)
			return 65; //18	0.052	65	602.201 Mb
		else if (numBalls <= 1 << 20)
			return 61; //20	0.06	61	2544.7 Mb
		else if (numBalls <= 1 << 22)
			return 68; //22	0.051	68	10748.8 Mb
		else if (numBalls <= 1 << 24)
			return 69; //24	0.051	69	45183 Mb
		
		return 69;
	}

	u64 SimpleIndex::getNumBins(u64 numBalls) {

		if (numBalls <= 1 << 8)
			return 11;
		else if (numBalls <= 1 << 10)
			return 56;
		else if (numBalls <= 1 << 12)
			return 204;
		else if (numBalls <= 1 << 14)
			return 868; //0.053	868	62	267514128 bits = 33.4393 Mb 
		else if (numBalls <= 1 << 16)
			return 3801; //0.058	3801	60	1134264012 bits = 141.783 Mb 
		else if (numBalls <= 1 << 18)
			return 13631; //0.052	13631	65	4817604330 bits = 602.201 Mb
		else if (numBalls <= 1 << 20)
			return 62914; //0.06	62914	61	20357586292 bits = 2544.7 Mb
		else if (numBalls <= 1 << 22)
			return 213909; //0.051	213909	68	85990134546 bits = 10748.8 Mb
		else if (numBalls <= 1 << 24)
			return 855638; //0.051	855638	69	361464273100 bits = 45183 Mb

		return 855638;
	}

    void SimpleIndex::init(u64 numBalls, u64 statSecParam)
    {
		mNumBins = getNumBins(numBalls);
		mMaxBinSize = getMaxBinSize(numBalls);

		mMaxBinSize += 1; //add 1 bot for PSU

		mBins.resize(mNumBins);
    }


    void SimpleIndex::insertItems(span<block> items, u64 numThreads)
    {

        AES hasher(mHashSeed);
		u64 inputSize = items.size();
		std::mutex mtx;
		const bool isMultiThreaded = numThreads > 1;

		auto routineHashing = [&](u64 t)
		{
			u64 inputStartIdx = inputSize * t / numThreads;
			u64 tempInputEndIdx = (inputSize * (t + 1) / numThreads);
			u64 inputEndIdx = std::min(tempInputEndIdx, inputSize);

			for (u64 i = inputStartIdx; i < inputEndIdx; ++i)
			{
				block temp = hasher.ecbEncBlock(items[i]);
				u64 addr = *(u64*)&temp % mNumBins;

					if (isMultiThreaded)
					{
						std::lock_guard<std::mutex> lock(mtx);
						mBins[addr].items.push_back(items[i]);
					}
					else
					{
						mBins[addr].items.push_back(items[i]);
					}

			}
		};

		std::vector<std::thread> thrds(numThreads);
		for (u64 i = 0; i < u64(numThreads); ++i)
		{
			thrds[i] = std::thread([=] {
				routineHashing(i);
			});
		}

		for (auto& thrd : thrds)
			thrd.join();
		
		//add a default block 
			for (u64 i = 0; i < mNumBins; ++i)
			{
			
				mBins[i].mBinRealSizes = mBins[i].items.size();
			}

    }

	
}
