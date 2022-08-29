import math


class SK42WGS:
    #Converts coordinats from SK42 to WGS system.
    #based on https://gis-lab.info/qa/wgs84-sk42-wgs84-formula.html
    math.pi = 3.14159265358979
    ro = 206264.8062
    aP = 6378245
    alP = 1 / 298.3
    e2P = 2 * alP - math.pow(alP, 2)
    aW = 6378137
    alW = 1 / 298.257223563
    e2W = 2 * alW - math.pow(alW, 2)
    a = float((aP + aW) / 2)
    e2 = float((e2P + e2W) / 2)
    da = aW - aP
    de2 = e2W - e2P
    dx = 23.92
    dy = -141.27
    dz = -80.9
    wx = 0
    wy = 0
    wz = 0
    ms = 0

    @staticmethod
    def WGS84_SK42_Lat(Bd, Ld, H):
        return Bd - SK42WGS.dB(Bd, Ld, H) / 3600

    @staticmethod
    def SK42_WGS84_Lat(Bd, Ld, H):
        return Bd + SK42WGS.dB(Bd, Ld, H) / 3600

    @staticmethod
    def WGS84_SK42_Long(Bd, Ld, H):
        return Ld - SK42WGS.dL(Bd, Ld, H) / 3600

    @staticmethod
    def SK42_WGS84_Long(Bd, Ld, H):
        return Ld + SK42WGS.dL(Bd, Ld, H) / 3600

    @staticmethod
    def dB(Bd, Ld, H):
        B = Bd * math.pi / 180
        L = Ld * math.pi / 180
        M = SK42WGS.a * (1 - SK42WGS.e2) / math.pow((1 - SK42WGS.e2 * math.pow(math.sin(B), 2)), 1.5)
        N = SK42WGS.a * math.pow((1 - SK42WGS.e2 * math.pow(math.sin(B), 2)), -0.5)
        return SK42WGS.ro / (M + H) * (N / SK42WGS.a * SK42WGS.e2 * math.sin(B) * math.cos(B) * SK42WGS.da
               + (math.pow(N, 2) / math.pow(SK42WGS.a, 2) + 1) * N * math.sin(B) * math.cos(B) * SK42WGS.de2 / 2
               - (SK42WGS.dx * math.cos(L) + SK42WGS.dy * math.sin(L)) * math.sin(B) + SK42WGS.dz * math.cos(B)) \
               - SK42WGS.wx * math.sin(L) * (1 + SK42WGS.e2 * math.cos(2 * B)) \
               + SK42WGS.wy * math.cos(L) * (1 + SK42WGS.e2 * math.cos(2 * B)) \
               - SK42WGS.ro * SK42WGS.ms * SK42WGS.e2 * math.sin(B) * math.cos(B)

    @staticmethod
    def dL(Bd, Ld, H):
        B = Bd * math.pi / 180
        L = Ld * math.pi / 180
        N = SK42WGS.a * math.pow((1 - math.pow(SK42WGS.e2 * math.sin(B), 2)), -0.5)
        return SK42WGS.ro / ((N + H) * math.cos(B)) * (-SK42WGS.dx * math.sin(L) + SK42WGS.dy * math.cos(L)) \
               + math.tan(B) * (1 - SK42WGS.e2) * (SK42WGS.wx * math.cos(L) + SK42WGS.wy * math.sin(L)) - SK42WGS.wz

    @staticmethod
    def WGS84Alt(Bd, Ld, H):
        B = Bd * math.pi / 180
        L = Ld * math.pi / 180
        N = SK42WGS.a * math.pow((1 - SK42WGS.e2 * math.pow(math.sin(B), 2)), -0.5)
        dH = -SK42WGS.a / N * SK42WGS.da + N * math.pow(math.sin(B), 2) * SK42WGS.de2 / 2 \
             + (SK42WGS.dx * math.cos(L) + SK42WGS.dy * math.sin(L)) * math.cos(B) + SK42WGS.dz * math.sin(B) \
             - N * SK42WGS.e2 * math.sin(B) * math.cos(B) \
             * (SK42WGS.wx / SK42WGS.ro * math.sin(L) - SK42WGS.wy / SK42WGS.ro * math.cos(L)) \
             + (math.pow(SK42WGS.a, 2) / N + H) * SK42WGS.ms
        return H + dH


print(SK42WGS.SK42_WGS84_Long(0, 0, 0), SK42WGS.SK42_WGS84_Lat(0, 0, 0), SK42WGS.WGS84Alt(0, 0, 0))