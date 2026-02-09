"""Loads and discovers FATES test classes"""

import importlib
import pkgutil
import tests.functional as functional_path
from framework.functional_test import FunctionalTest


def discover_tests():
    """Scans the tests/functional directory and imports everything"""
    
    registry = {}
    # walk through all subdirectories in tests/functional
    for _, module_name, _ in pkgutil.walk_packages(
        functional_path.__path__, functional_path.__name__ + "."
    ):
        module = importlib.import_module(module_name)

        # look for any class inside the module that is a subclass of FunctionalTest
        for _, obj in vars(module).items():
            if (
                isinstance(obj, type)
                and issubclass(obj, FunctionalTest)
                and obj is not FunctionalTest
            ):
                # use the 'name' attribute defined in your class to register it
                registry[obj.name] = obj

    return registry


def get_test_instances(config_dict: dict) -> dict:
    """
    Transforms the raw config dictionary into a dictionary of
    live test objects.
    """
    # trigger test discovery
    discover_tests()

    # map subclasses
    available_classes = {}
    for base in [FunctionalTest]:
        for cls in base.__subclasses__():
            # ensure we only pick concrete classes, not the base classes themselves
            if hasattr(cls, "name"):
                available_classes[cls.name] = cls

    # instatiate test classes
    test_instances = {}
    for name, attributes in config_dict.items():
        if name in available_classes:
            # pass attributes dict to the constructor
            test_instances[name] = available_classes[name](name, attributes)
        else:
            print(
                f"WARNING: No Python class found for test '{name}'. "
                f"Check that 'name = \"{name}\"' is defined in your test class."
            )

    return test_instances

def validate_test_configs(test_instances: dict):
    """Checks that all external dependencies (like DATM files) exist before running any
        tests.

    Args:
        test_instances (dict): dictionary of test class instances

    Raises:
        FileNotFoundError: Can't find file
    """
    
    missing_assets = []
    for name, test in test_instances.items():
        # only check tests that actually have a driver file defined
        if test.datm_file:
            if not test.datm_file.exists():
                missing_assets.append(f"[{name}] Driver missing: {test.datm_file}")

    if missing_assets:
        error_msg = "\n".join(missing_assets)
        raise FileNotFoundError(f"Pre-run validation failed:\n{error_msg}")
    