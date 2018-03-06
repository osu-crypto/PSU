#include "stdafx.h"
#ifdef  _MSC_VER
#include "CppUnitTest.h"
#include "PSU_Tests.h"
#include "NcoOT_Tests.h"
#include "Common.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace tests_libOTe
{
    TEST_CLASS(nOPRF_Tests)
    {
    public:

       
        TEST_METHOD(Poly_IKNP_TestVS)
        {
            InitDebugPrinting();
			Poly_IKNP_Test_Impl();
        }
		 

		TEST_METHOD(PMT_TestVS)
		{
			InitDebugPrinting();
			PMT_Test_Impl();
		}

    };
}
#endif