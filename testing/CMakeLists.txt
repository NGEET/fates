# This is where you add specific test directories

## Functional tests
add_subdirectory(tests/functional/allometry fates_allom_ftest)
add_subdirectory(tests/functional/fire/fuel fates_fuel_ftest)
add_subdirectory(tests/functional/fire/ros fates_ros_ftest)
add_subdirectory(tests/functional/patch fates_patch_ftest)
add_subdirectory(tests/functional/fire/mortality fates_firemort_ftest)

## Unit tests
add_subdirectory(tests/unit/fire_weather_test fates_fire_weather_utest)
add_subdirectory(tests/unit/fire_fuel_test fates_fire_fuel_utest)
add_subdirectory(tests/unit/sort_cohorts_test fates_sort_cohorts_utest)
add_subdirectory(tests/unit/insert_cohort_test fates_insert_cohort_utest)
add_subdirectory(tests/unit/validate_cohorts_test fates_validate_cohorts_utest)
add_subdirectory(tests/unit/count_cohorts_test fates_count_cohorts_utest)
add_subdirectory(tests/unit/fire_equations_test fates_fire_equations_utest)
add_subdirectory(tests/unit/quadratic_roots_test fates_quadratic_roots_utest)
add_subdirectory(tests/unit/great_circle_test fates_great_circle_utest)
