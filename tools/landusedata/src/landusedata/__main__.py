from landusedata._main import main

# Gaurd against import time side effects
if __name__ == '__main__':
    raise SystemExit(main())
