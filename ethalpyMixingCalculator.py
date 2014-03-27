import numpy as np
import matplotlib.pylab as plt

class hmix():
    def __init__(self):
        self.A_phiStar = 0.0
        self.A_nws13 = 0.0
        self.A_Vm23 = 0.0
        self.A_name = ''
        self.B_phiStar = 0.0
        self.B_nws13 = 0.0
        self.B_Vm23 = 0.0
        self.B_name = ''
        self.deH_A_partial_infDilute = 0.0
        self.P = 0.0
        self.RP = 0.0
        self.QP = 9.4
        self.e = 1.0
        self.elementName = {}
        self.elementPhiStar = {}
        self.elementNWS13 = {}
        self.elementVM23 = {}
        self.elementRP = {}
        self.elementTRAN = {}
        self.Avogardro = 6.02E23 # unit /mole
        self.xA = np.linspace(0.001,0.999,200)
        self.xB = 1.0 - self.xA        
        self.__printHead()        
        self.__readDatabase()
        
    def __printHead(self):
        print '-------------------------------------------------'
        print '-Calculate heat of mixing of binary metal system-'
        print '-            based on Miedema scheme            -'
        print '-              Writen by Jun Ou                 -'
        print '-------------------------------------------------'
        return

    def __readDatabase(self):
        try:
            fReadDatabase = open('./database.dat','r')
            for line in fReadDatabase:
                elementTmp = line.replace('\n',' ').replace(' ','').split(',')
                
                #print elementTmp
                nameTmp = elementTmp[0]
                self.elementName[nameTmp] = elementTmp[0]                
                self.elementPhiStar[nameTmp] = elementTmp[1]                
                self.elementNWS13[nameTmp] = elementTmp[2]
                self.elementVM23[nameTmp] = elementTmp[3]
                self.elementRP[nameTmp] = elementTmp[4]
                self.elementTRAN[nameTmp] = elementTmp[5]
            print 'succesful initialization of the database.'    
        except ValueError:
            print 'database Error'
        return
    
    # get input for the two components    
    def inputElement(self): 
        self.A_name = raw_input('Please input the first component: ')
        self.B_name = raw_input('Please input the second component: ')
        #self.A_name = 'Ti'
        #self.B_name = 'Al'
        return
    
    def calRP(self):
        if self.elementTRAN[self.A_name]=='T' and self.elementTRAN[self.B_name]=='T':
            self.RP = 0.0
        elif self.elementTRAN[self.A_name]=='N' and self.elementTRAN[self.B_name]=='N':
            self.RP = 0.0
        else:
            self.RP = float(self.elementRP[self.A_name])*float(self.elementRP[self.B_name])*0.73
        return
    
    def assiginP(self):
        if self.elementTRAN[self.A_name]=='T' and self.elementTRAN[self.B_name]=='T':
            self.P = 0.147
        elif self.elementTRAN[self.A_name]=='N' and self.elementTRAN[self.B_name]=='N':
            self.P = 0.111
        else:
            self.P = 0.128            
        return    
    
    def calHmix(self):
        self.A_phiStar = float(self.elementPhiStar[self.A_name])
        self.B_phiStar = float(self.elementPhiStar[self.B_name])
        self.A_nws13 = float(self.elementNWS13[self.A_name])
        self.B_nws13 = float(self.elementNWS13[self.B_name])
        self.A_Vm23 = float(self.elementVM23[self.A_name])
        self.B_Vm23 = float(self.elementVM23[self.B_name])
        dePhi = self.A_phiStar-self.B_phiStar
        deNws13 = self.A_nws13-self.B_nws13
        self.xAs = self.xA*self.A_Vm23/(self.xA*self.A_Vm23+self.xB*self.B_Vm23)
        self.fxs = self.xAs*(1.0-self.xAs)
        self.g = 2.0*(self.xA*self.A_Vm23+self.xB*self.B_Vm23)/(1.0/self.A_nws13+1.0/self.B_nws13)
        self.deHmix = self.Avogardro*self.fxs*self.g*self.P*(-self.e*(dePhi)**2+self.QP*(deNws13)**2-self.RP)*1.60217657E-22
        self.deH_A_partial_infDilute =2.0*self.A_Vm23/(1.0/self.A_nws13+1.0/self.B_nws13)*self.Avogardro*self.P*(-self.e*(dePhi)**2+self.QP*(deNws13)**2-self.RP)*1.60217657E-22
        print 'Calculation done.'
        return
    
    def report(self):
        print ''
        print '-------------------------------------------------'
        print '-------------------- report ---------------------'
        print 'Two components: ', self.A_name, '(',self.elementTRAN[self.A_name], ') and ', self.B_name, '(',self.elementTRAN[self.B_name],')'
        print 'Phi of', self.A_name, self.A_phiStar 
        print 'Phi of', self.B_name, self.B_phiStar
        print 'n(ws)^1/3 of', self.A_name, self.A_nws13
        print 'n(ws)^1/3 of', self.B_name, self.B_nws13
        print 'Vm^2/3 of', self.A_name, self.A_Vm23
        print 'Vm^2/3 of', self.B_name, self.B_Vm23
        
        print 'P:', self.P
        print 'R/P:', self.RP
        print 'Q0/P:',self.QP
        print 'deHA at infinite dilution', self.deH_A_partial_infDilute
        
        print '-------------------------------------------------'
        return
              
    def plot(self):
        fooPlot = raw_input('Plot a figure? y (default)/n: ')
        if fooPlot=='n' or fooPlot=='N':
            print 'No is chosen. No figure is produced.'       
        else:
            print 'A figure has been produced.'
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('Mole Fraction of '+str(self.A_name)+' in '+str(self.A_name)+'-'+str(self.B_name))
            ax.set_ylabel(r'$\Delta H_{mix}$ [kJ/mole]')
            ax.set_xlim(0,1)
            ax.plot(self.xA,self.deHmix)
            plt.show()
            plt.close(fig) 
        return
              
def main():
    H = hmix()
    while True:
        H.inputElement()
        H.calRP()
        H.assiginP()
        H.calHmix()
        H.report()
        H.plot()
    return
    
if __name__=='__main__':
    main()
