from ctypes import *
import numpy as np
import time
mydll = cdll.LoadLibrary(r".\test.dll")
z = range(1024)
x = np.array(z)
out = np.zeros([1024])
delay_buf = np.zeros([1024])
x_len = len(x)
prev_i =c_double(1.1)
print(delay_buf)
mydll.test(c_void_p(delay_buf.ctypes.data))
mydll.test(c_void_p(delay_buf.ctypes.data))

for i in range(10):
    t0 = time.clock()
    mydll.chorus(c_void_p(x.ctypes.data),c_void_p(out.ctypes.data),\
                c_void_p(delay_buf.ctypes.data),x_len,addressof(prev_i))
    print( time.clock() - t0)
encoded =  out.astype('float32').tobytes()
print(delay_buf[:20])
print('asdf:',prev_i)
