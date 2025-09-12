import sys
import os

def _run():
    # Simple test runner for environments without pytest
    failures = 0
    # import tests directly
    try:
        from tests import test_simulator_basic as t
    except Exception as e:
        print('ERROR: failed to import tests:', e)
        return 2
    for name in dir(t):
        if name.startswith('test_'):
            fn = getattr(t, name)
            print(f'RUN {name}...', end='')
            try:
                fn()
                print(' OK')
            except Exception as e:
                print(' FAIL')
                import traceback
                traceback.print_exc()
                failures += 1
    if failures:
        print(f'{failures} test(s) failed')
        return 1
    print('ALL TESTS PASSED')
    return 0

if __name__ == '__main__':
    sys.exit(_run())
