add_test( MinkowskiTests.WorkingPoints /Users/jacklewis/Documents/work/year3/DMet/build/test/Google_tests_run [==[--gtest_filter=MinkowskiTests.WorkingPoints]==] --gtest_also_run_disabled_tests)
set_tests_properties( MinkowskiTests.WorkingPoints PROPERTIES WORKING_DIRECTORY /Users/jacklewis/Documents/work/year3/DMet/build/test SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set( Google_tests_run_TESTS MinkowskiTests.WorkingPoints)
