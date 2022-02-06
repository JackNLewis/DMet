add_test( MinkowskiTests.WorkingInRange /Users/jacklewis/Documents/work/year3/DMet/build/test/Google_tests_run [==[--gtest_filter=MinkowskiTests.WorkingInRange]==] --gtest_also_run_disabled_tests)
set_tests_properties( MinkowskiTests.WorkingInRange PROPERTIES WORKING_DIRECTORY /Users/jacklewis/Documents/work/year3/DMet/build/test SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
add_test( MinkowskiTests.SamePoints /Users/jacklewis/Documents/work/year3/DMet/build/test/Google_tests_run [==[--gtest_filter=MinkowskiTests.SamePoints]==] --gtest_also_run_disabled_tests)
set_tests_properties( MinkowskiTests.SamePoints PROPERTIES WORKING_DIRECTORY /Users/jacklewis/Documents/work/year3/DMet/build/test SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
add_test( ChebyshevTests.WorkingPoints /Users/jacklewis/Documents/work/year3/DMet/build/test/Google_tests_run [==[--gtest_filter=ChebyshevTests.WorkingPoints]==] --gtest_also_run_disabled_tests)
set_tests_properties( ChebyshevTests.WorkingPoints PROPERTIES WORKING_DIRECTORY /Users/jacklewis/Documents/work/year3/DMet/build/test SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set( Google_tests_run_TESTS MinkowskiTests.WorkingInRange MinkowskiTests.SamePoints ChebyshevTests.WorkingPoints)
