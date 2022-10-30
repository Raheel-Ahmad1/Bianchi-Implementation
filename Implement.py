from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math 
import random


class Bianchi:
    # Calculates the delay and jitter of IEEE 802.11 WLAN basic scheme according to:
    # G.Bianchi, "Performance analysis of the IEEE 802.11 distributed coordination function," in IEEE
    # Journal on Selected Areas in Communications, vol. 18, no. 3, pp. 535 - 547, March 2000.
    # doi: 10.1109 / 49.840210
    def __init__(self, bitrate, n, ACK, SIFS, slot, DIFS, E_P, E_P_star, W, m, H, prop_delay):
        # INPUT:
        # bitrate: raw bitrate in bps
        # n: number of STAs
        # ACK: ACK length in bits
        # SIFS: SIFS duration in seconds
        # slot: slot duration in seconds
        # DIFS: DIFS duration in seconds
        # E_P: average packet payload size in bits
        # E_P_star:  average  length  of  the  longest packet payload involved in a collision in bit (for an example, see eq. 16 in the paper)
        # W: Minimum contention window size in slots
        # m: Retry limit
        # H: Header size in bits
        # prop_delay: Propagation delay in seconds
        #
        # OUTPUT:
        #
        # p: the probability of a collision seen by a packet being transmitted on the channel
        # t: the probability that a station transmits in a randomly chosen slot time
        # Ps: the probability that a transmission occurring on the channel is successful is given by the probability that exactly one station transmits on the channel, conditioned on the fact that at least one station transmits
        # P_tr: the probability that there is at least one transmission in the considered slot time
        # T_s: the average time the channel is sensed busy, in seconds
        # T_c: the average time the channel is sensed busy by each station during a collision in seconds.

        # independent varz
        self.bitrate = 56E6
        self.n = 30

        self.b00 = 0

        self.ACK = 14 * 8
        self.SIFS = 10E-6
        self.slot = 9E-6
        self.DIFS = 3 * 9E-6
        self.H = 2 * 8

        self.E_P = 80
        self.E_P_star = 100

        self.W = 12
        self.m = 6

        self.prop_delay = 0

        # dependent varz
        self.p = 0
        self.t = 0
        self.Ps = 0
        self.P_tr = 0
        self.T_s = 0
        self.T_c = 0
        self.S = 0
        self.E_s=0
        self.P_e=0
        self.P_c=0
        
        self.bitrate = bitrate
        self.n = n
        self.ACK = ACK
        self.SIFS = SIFS
        self.slot = slot
        self.DIFS = DIFS
        self.E_P = E_P
        self.E_P_star = E_P_star
        self.W = W
        self.m = m
        self.H = H
        self.prop_delay = prop_delay

        self.calculate_p_t()
        self.calculate_Ptr()
        self.calculate_Ps()
        self.calculate_Ts()
        self.calculate_Tc()
        self.calculate_S()
        self.calculate_b00()
        self.calculate_Es()

    def calculate_b00(self):
        self.b00 = 2.0 * (1.0 - 2.0 * self.p) * (1.0 - self.p) / ((1.0 - 2.0 * self.p) * (self.W + 1) + self.p * self.W * (1.0 - (2.0 * self.p)**self.m))

    def calculate_b(self, i, k):
        if i < self.m:
            W_i = 2.0 ** i * self.W
            b_i0 = self.p ** i * self.b00
        else:
            W_i = 2.0 ** self.m * self.W
            b_i0 = self.p ** self.m / (1 - self.p) * self.b00

        return (W_i - k) / W_i * b_i0

    def check_p_t(self):
        c1 = self.p - 1.0 + (1.0 - self.t) ** (self.n - 1.0) <= 1.49012e-08
        my_sum = 0.0
        for i in range(0, self.m):
            my_sum += (2.0 * self.p) ** i
        c2 = 2.0 / (1.0 + self.W + self.p * self.W * my_sum) - self.t <= 1.49012e-08
        return c1 and c2

    def calculate_p_t(self):
        def equations(x):
            p, t = x
            my_sum = 0.0
            for i in range(0, self.m):
                my_sum += (2.0 * p) ** i

            return ( p - 1.0 + (1.0 - t) ** (self.n - 1.0), 2.0 / (1.0 + self.W + p * self.W * my_sum) - t )

        self.p, self.t = fsolve(equations, (0.1, 0.1))

        print ("System solved, error (p, tau): " + str(equations((self.p, self.t))))

    def calculate_Ps(self):
        self.Ps = self.n * self.t * (1.0 - self.t) ** (self.n - 1) / self.P_tr

    def calculate_Ptr(self):
        self.P_tr = 1.0 - (1.0 - self.t) ** self.n

    def calculate_Ts(self):
        self.T_s = self.H / self.bitrate + self.E_P / self.bitrate + self.SIFS + self.prop_delay + self.ACK / self.bitrate + self.DIFS + self.prop_delay

    def calculate_Tc(self):
        self.T_c = self.H / self.bitrate + self.E_P_star / self.bitrate + self.DIFS + self.prop_delay
        
    def calculate_Es(self): ##Average slot duration
        self.P_e=1-self.P_tr
       # self.P_c=self.P_tr-self.Ps
        self.E_s = self.P_e*9E-6 + self.Ps*self.T_s*self.P_tr + self.P_tr*self.T_c*(1-self.Ps)
        P=[]
        Dly=[]
        Jit=[]
        E_D=0
        E_D2=0
        for i in range(0, 7):
            E_U0=self.T_c+self.E_s*((self.W-1)/2)
            w_i=0
            for k in range(0,i):
                w_i=w_i+((2.0 ** k * self.W)-1)/2
                
            E_Uj=(i+1)*self.T_c+self.E_s*w_i
            Prob=((1-self.p)*self.p**i)/((1-self.p**(self.m+1)))
            E_D1=0
            
            for p in range(0,(2 ** i * self.W)-1):
                E_D1=E_D1+(self.T_s+p*self.E_s+E_Uj+E_U0)**2
        
            P.append(Prob)
            print('Probability',Prob)
            a=self.calculate_b(i,self.m)
            delay=self.T_s+a*self.E_s+E_Uj+E_U0
            Dly.append(delay*10E3)
            E_D=E_D+(self.T_s+i*self.T_c+self.E_s*w_i)*(((1-self.p)*self.p**i)/(1-self.p**(self.m+1)))
            E_D2=E_D2+((1-self.p)*self.p**i)/((1-self.p**(self.m+1))*2.0 ** i * self.W)*E_D1
            
           # Jit.append(jitter)
        jitter=E_D2-E_D**2   
        print('Delay',Dly)
        print('jitter==',math.sqrt(jitter))
        random_number = random.choices(Dly, P)

       # print('random',random_number[0])
        a=[0,1,2,3,4,5,6]
       
      
        plt.plot(a,P)
        plt.xlabel("Stage ")
        plt.ylabel("Probability")
        plt.figure()
        plt.plot(Dly,P)
        
        plt.xlabel("End-to-End Packet Delay")
        plt.ylabel("Probability")
        plt.tight_layout()
        plt.show()
   
        file = open("myfile.txt", "w")
        val=Dly[0]
        for i in range(0,300000):
           # random_number=random.choices(Dly, P)
           # random_number=( "{0:.3f}".format(random_number[0]/1000))
           # file.write(str(random_number))
           # file.write(" ")
            val=val+0.001
            word=( "{0:.6f}".format(val/1000))
            file.write(str(word))
            file.write(" ")
        file.close()
    def calculate_S(self):
        self.S = self.Ps * self.P_tr * (self.E_P / self.bitrate) / (
        (1.0 - self.P_tr) * self.slot + self.P_tr * self.Ps * self.T_s + self.P_tr * (1.0 - self.Ps) * self.T_c)


class uniform_helper:
    # Calculates the parameters for a uniformly distributed packet length size between P_min and P_max (in bits)

    # independent varz
    P_min = 0
    P_max = 0

    #dependent varz
    E_P = 0
    E_P_star = 0

    def __init__(self, P_min, P_max):
        self.P_max = P_max
        self.P_min = P_min
        self.calculate_EP()
        self.calculate_E_P_star()

    def calculate_EP(self):
        self.E_P = (self.P_max - self.P_min) / 2.0

    def calculate_E_P_star(self):
        self.E_P_star = ( (self.P_max - self.P_min) ** 2.0 - 1 ) / (self.P_max - self.P_min)


if __name__ == "__main__":
    # Validation: see Fig.6 Bianchi paper

    bitrate = 56E6
    ACK = 112 + 128
    SIFS = 10E-6
    slot = 9E-6
    DIFS = 27E-6
    E_P = 8184
    E_P_star = E_P
    WW = [32, 128]
    mm = [3, 5]
    H = 272 + 128
    prop_delay = 0

    fig = plt.figure()

    ax = fig.gca()

    nn = range(5, 50)
    B = Bianchi(bitrate, 25, ACK, SIFS, slot, DIFS, E_P, E_P_star, 32, 5, H, prop_delay)
