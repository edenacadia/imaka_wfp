## CorrClock
## Eden McEwen
## October 2021
# Timer class for correlating imaka AOCB

import time

class CorrClock(object):
    def __init__(self, acor_t = 0, xcor_t = 0):
        self.acor_init = 0
        self.xcor_init = 0
        self.acor_end = acor_t
        self.xcor_end = xcor_t
        self.acor_par = False
        self.xcor_par = False
        
    def acor_start(self, par=False):
        self.acor_par = par
        self.acor_init = time.time()
        
    def acor_stop(self):
        self.acor_end = time.time()
        
    def acor_time(self):
        return self.acor_end - self.acor_init
        
    def xcor_start(self, par=False):
        self.xcor_par = par
        self.xcor_init = time.time()
        
    def xcor_stop(self):
        self.xcor_end = time.time()
        
    def xcor_time(self):
        return self.xcor_end - self.xcor_init