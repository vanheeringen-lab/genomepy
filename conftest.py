import py.test

def pytest_addoption(parser):
    parser.addoption('--slow', action='store_true', default=False,
                      help='Also run slow tests')

def pytest_runtest_setup(item):
    """Skip tests if they are marked as slow and --slow is not given"""
    if getattr(item.obj, 'slow', None) and not item.config.getvalue('slow'):
        py.test.skip('slow tests not requested')
