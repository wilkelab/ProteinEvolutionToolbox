import re

class DSSPData:
  def __init__(self):
      self.num    = []
      self.resnum = []
      self.moltyp = []
      self.aa     = []
      self.struct = []
      self.bp1    = []
      self.bp2    = []
      self.acc    = []
      self.h_nho1 = []
      self.h_ohn1 = []
      self.h_nho2 = []
      self.h_ohn2 = []
      self.tco    = []
      self.kappa  = []
      self.alpha  = []
      self.phi    = []
      self.psi    = []
      self.xca    = []
      self.yca    = []
      self.zca    = []

  def parseDSSP(self, file):
    input_handle = open(file, 'r')

    line_num = 0
    start=False
    for line in input_handle:
  
      if( re.search('#', line) ):
        start=True
        continue

      if( start ):
        self.num.append(    line[0:5].strip() )
        self.resnum.append( line[5:10].strip() )
        self.moltyp.append( line[10:12].strip() )
        self.aa.append(     line[12:14].strip() )
        self.struct.append( line[14:25] )
        self.bp1.append(    line[25:29].strip() )
        self.bp2.append(    line[29:34].strip() )
        self.acc.append(    line[34:38].strip() )
        self.h_nho1.append( line[38:50].strip() )
        self.h_ohn1.append( line[50:61].strip() )
        self.h_nho2.append( line[61:72].strip() )
        self.h_ohn2.append( line[72:83].strip() )
        self.tco.append(    line[83:91].strip() )
        self.kappa.append(  line[91:97].strip() )
        self.alpha.append(  line[97:103].strip() )
        self.phi.append(    line[103:109].strip() )
        self.psi.append(    line[109:115].strip() )
        self.xca.append(    line[115:122].strip() )
        self.yca.append(    line[122:129].strip() )
        self.zca.append(    line[129:136].strip() )

  def getResnums(self):
    return self.resnum
  def getChainType(self):
    return self.moltyp
  def getAAs(self):
    return self.aa
  def getSecStruc(self):
    return self.struct
  def getBP1(self):
    return self.bp1
  def getBP2(self):
    return self.bp2
  def getACC(self):
    return self.acc
  def getH_NHO1(self):
    return self.h_nho1
  def getH_NHO2(self):
    return self.h_nho2
  def getH_OHN1(self):
    return self.h_ohn1
  def getH_OHN2(self):
    return self.h_ohn2
  def getTCO(self):
    return self.tco
  def getKAPPA(self):
    return self.kappa
  def getALPHA(self):
    return self.alpha
  def getPHI(self):
    return self.phi
  def getPSI(self):
    return self.psi
  def getX(self):
    return self.xca
  def getY(self):
    return self.yca
  def getZ(self):
    return self.zca



