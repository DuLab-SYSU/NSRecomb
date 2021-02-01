import functools
import time


def execute_time(func):  # 定义装饰器
    # 此时调用f()时，实际上调用的是wrapper函数，返回值一样，但名字变了。
    @functools.wraps(func)  # 把名字改回来
    def wrapper(*args, **kw):
        print('Start Task ======================>')
        start_time = time.time()
        res = func(*args, **kw)
        end_time = time.time()
        print("===========> The time is %10.4s s." % (end_time-start_time))
        return res
    return wrapper
