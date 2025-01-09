# This is where you add specific test directories

## Functional tests
add_subdirectory(functional_testing/allometry fates_allom_ftest)
add_subdirectory(functional_testing/math_utils fates_math_ftest)
add_subdirectory(functional_testing/fire/fuel fates_fuel_ftest)
add_subdirectory(functional_testing/fire/ros fates_ros_ftest)

## Unit tests
add_subdirectory(unit_testing/fire_weather_test fates_fire_weather_utest)
add_subdirectory(unit_testing/fire_fuel_test fates_fire_fuel_utest)
