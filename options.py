from math import exp, sqrt
from datetime import datetime, timedelta, date
from numpy import is_busday, busday_count

# Helpers

def discount(amt, days, rfr):
    # positive days means future payout (so value goes down)
    return amt * exp(-rfr*(days/365.))

def exercise_call(stock,strike):
    return stock - strike

def exercise_put(stock, strike):
    return strike - stock

# Option class. Parent to BinOption and TrinOption
class Option:
    def __init__(self,info):
        self.today = datetime.strptime(info['today'],'%m/%d/%Y')
        self.type = info['type'] # 'call' or 'put'
        self.price = info['price']
        self.strike = info['strike']
        self.vola = info['vola']
        self.expi = datetime.strptime(info['expi'],'%m/%d/%Y')
        self.rfr = info['rfr']
        try:
            # dividends should be in form (date, amt)
            self.divs = [(datetime.strptime(date,'%m/%d/%Y'),amt) for (date, amt) in info['divs']]
            self.divs = [div for div in self.divs if div[0]>self.today]
            self.divs.sort(key=lambda x: x[0])
        except KeyError:
            self.divs = []
        try:
            # input levels
            self.periods = info['periods']
        except KeyError:
            self.periods = 3
        self.tstep = (self.expi-self.today).days / float(self.periods)

# Helpers for binomial tree option pricer

# calculate high price multipliers (low_mult = 1/high_mult)
def high_mult(vola,tstep):
    return exp(vola * sqrt(tstep/365.))

# calculate high price probability (low_p = 1-high_p)
def high_p(rfr, vola, tstep):
    u = high_mult(vola,tstep)
    d = 1./u
    return (exp(rfr*tstep/365.) - d)/(u-d)

# BinOption class. May make it a subclass of an Option
class BinOption(Option):
    # generates binomial price tree
    def pricetree(self):
        if self.type == 'call':
            exercise = exercise_call
        elif self.type == 'put':
            exercise = exercise_put
        up = high_mult(self.vola,self.tstep)
        p_up = high_p(self.rfr,self.vola,self.tstep)
        p_down = 1-p_up
        # pricetree holds option value at each node
        pricetree = [0]*((self.periods+1)*(self.periods+2)/2)
        # calculate initial price, minus initial PV of divs
        start_div_val = 0
        for div in self.divs:
            start_div_val += discount(div[1],(div[0]-self.today).days,self.rfr)
        adj_price = self.price - start_div_val
        # generate pricetree
        day = self.expi
        cur_div_val = 0
        last_div = len(self.divs)-1
        for step in xrange(self.periods,-1,-1):
            # find exercise values via ups/downs, add back PV of divs
            firstNode = step*(step+1)/2
            childNode = (step+1)*(step+2)/2
            for i in range(step+1):
                # exercise value
                price = adj_price * up**(2*i-step)
                price += cur_div_val
                exer_val = exercise(price, self.strike)
                # binomial value
                bin_val = 0
                if step < self.periods:
                    weighted_avg = p_down*pricetree[childNode+i]+p_up*pricetree[childNode+i+1]
                    bin_val = discount(weighted_avg,self.tstep,self.rfr)
                # final value
                pricetree[firstNode+i] = max(bin_val,exer_val,0)
            # recalculate PV of future dividends at this step
            day -= timedelta(days=self.tstep)
            if (last_div>=0) and (self.divs[last_div][0]>day):
                cur_div_val+=self.divs[last_div][1]
                last_div -= 1
            cur_div_val = discount(cur_div_val,self.tstep,self.rfr)
        return pricetree

    def value(self):
        return self.pricetree()[0]

# Helpers for trinomial tree option pricer

def tri_high_mult(vola,tstep):
    return exp(vola*sqrt(2*tstep/365.))

def tri_low_mult(vola,tstep):
    return exp(-vola*sqrt(2*tstep/365.))

def tri_high_p(rfr,vola,tstep):
    const = vola*sqrt(tstep/365./2.)
    num = exp(rfr*tstep/365./2.) - exp(-const)
    den = exp(const) - exp(-const)
    return (num/den)**2

def tri_low_p(rfr,vola,tstep):
    const = vola*sqrt(tstep/365./2.)
    num = exp(const) - exp(rfr*tstep/365./2.)
    den = exp(const) - exp(-const)
    return (num/den)**2

class TrinOption(Option):
    # generates trinomial tree
    def pricetree(self):
        if self.type == 'call':
            exercise = exercise_call
        elif self.type == 'put':
            exercise = exercise_put
        up = tri_high_mult(self.vola,self.tstep)
        p_up = tri_high_p(self.rfr,self.vola,self.tstep)
        p_down = tri_low_p(self.rfr,self.vola,self.tstep)
        p_mid = 1 - p_up - p_down
        # pricetree holds option value at each node
        pricetree = [0]*(self.periods+1)**2
        # calculate initial price, minus initial PV of divs
        start_div_val = 0
        for div in self.divs:
            start_div_val += discount(div[1],(div[0]-self.today).days,self.rfr)
        adj_price = self.price - start_div_val
        # generate pricetree
        day = self.expi
        cur_div_val = 0
        last_div = len(self.divs)-1
        for step in xrange(self.periods,-1,-1):
            # find exercise values via ups/downs, add back PV of divs
            firstNode = step ** 2
            childNode = (step+1) ** 2
            for i in range(2*step+1):
                # exercise value
                price = adj_price * up**(i-step)
                price += cur_div_val
                exer_val = exercise(price, self.strike)
                # trinomial value
                trin_val = 0
                if step < self.periods:
                    weighted_avg = p_down*pricetree[childNode+i]+p_mid*pricetree[childNode+i+1]+p_up*pricetree[childNode+i+2]
                    trin_val = discount(weighted_avg,self.tstep,self.rfr)
                # final value
                pricetree[firstNode+i] = max(trin_val,exer_val,0)
            # recalculate PV of future dividends at this step
            day -= timedelta(days=self.tstep)
            if (last_div>=0) and (self.divs[last_div][0]>day):
                cur_div_val+=self.divs[last_div][1]
                last_div -= 1
            cur_div_val = discount(cur_div_val,self.tstep,self.rfr)
        return pricetree

    def value(self):
        return self.pricetree()[0]

# Example data
# divs = [('4/19/2018',2),('4/21/2018',2)]
# testInfo = {'divs':divs, 'today':'1/19/2018','type':'call','price':50,
#             'strike':50,'vola':0.4,'expi':'5/19/2018','rfr':0.09,'periods':200}
# testOp = BinOption(testInfo)
# print testOp.value()

