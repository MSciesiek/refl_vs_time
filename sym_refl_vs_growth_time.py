import numpy as np
import matplotlib.pyplot as plt


# %% funkcje
def Pmatrix(n, d, lamb):
    '''
    Macierz przejscia fali o dl 'lamb' przed material o grubosci 'd' i wsp.
    zalamania 'n'
    '''
    return np.matrix([[np.exp(-1j*2*np.pi/lamb*n*d), 0],
                      [0, np.exp(1j*2*np.pi/lamb*n*d)]])


def Dmatrix(n1, n2):
    '''
    Macierz przejscia przez granice osrodkow o wsp. zalamania 'n1' i 'n2' z
    osrodka n2 do n1
    '''
    t = 2*n1/(n1+n2)
    r = (n1-n2)/(n1+n2)
    return np.matrix([[1/t, r/t], [r/t, 1/t]])


def PD(t):
    '''
    Macierz przejscia dla odcinka wyhodowanego w t-tej sekundzie czasu
    '''
    if t > 0:
        if material[t] == 0:  # n_lo
            if material[t-1] == 0:
                return Pmatrix(n_lo, v_lo, wl)
            else:
                return Pmatrix(n_lo, v_lo, wl)*Dmatrix(n_lo, n_hi)
        elif material[t] == 1:  # n_hi
            if material[t-1] == 1:
                return Pmatrix(n_hi, v_hi, wl)
            else:
                return Pmatrix(n_hi, v_hi, wl)*Dmatrix(n_hi, n_lo)
        elif material[t] == 2:  # powietrze
            return Dmatrix(n_air, n_hi)
    else:  # t== 0, podloze
        return Dmatrix(n_hi, n_subst)


def M(t, _cache={}):
    '''
    Macierz przejscia dla struktury wyhodowanej do t-tej sekundy czasu
    wlacznie
    '''
    if t in _cache:
        return _cache[t]
    elif t == 0:
        return PD(0)
    else:
        return _cache.setdefault(t, PD(t)*M(t-1))

# %% parametry
# Parametry DBR
lamb0 = 800  # nm
DBR_bot = 26
DBR_mid = 16
DBR_top = 25

# grubosci odczytane z obrazka 0707_3_cd.jpg
d_hi = 85.4  # nm
d_lo = 97.9  # nm

n_air = 1
absorb = 0  # 1j*0.001  # <------------- !!!!!
n_hi = lamb0/(4*d_hi) + absorb
n_lo = lamb0/(4*d_lo)
n_subst = 3

# czasy odczytane z zeszytu
t_buffer = 3600  # s
t_lo = 172  # s
t_hi = 145  # s
t_cavity = 4*t_hi  # s

v_hi = d_hi/t_hi
v_lo = d_lo/t_lo
# badana dlugosc fali
wl = 800


# %% 'profil' n_lo, n_hi
material = []

test = 0
if test:
    mat0 = [0 for i in range(1000)]
    mat1 = [1 for i in range(1000)]  # interferencja w cienkiej warstwie
    mat2 = [0 for i in range(1000)]
    material = mat0 + mat1 + mat2
else:
    # 1 = 'hi', 0 = 'lo'
    for i in range(t_buffer):
        material.append(1)

    for i in range(DBR_bot-1):
        for t in range(t_lo):
            material.append(0)
        for t in range(t_hi):
            material.append(1)

    for t in range(t_lo):
        material.append(0)

    for t in range(t_cavity):
        material.append(1)

    for i in range(DBR_mid-1):
        for t in range(t_lo):
            material.append(0)
        for t in range(t_hi):
            material.append(1)

    for t in range(t_lo):
        material.append(0)

    for t in range(t_cavity):
        material.append(1)

    for i in range(DBR_top):
        for t in range(t_lo):
            material.append(0)
        for t in range(t_hi):
            material.append(1)

t_max = len(material)


# %% obliczenia
Ref = []

for t in range(t_max):
    if material[t] == 0:
        Mat = Dmatrix(n_air, n_lo)*M(t)
    else:
        Mat = Dmatrix(n_air, n_hi)*M(t)
    Ref.append((abs(Mat[1, 0])/abs(Mat[0, 0]))**2)

# outfile = "sym_DBR_refl_vs_delta_n.csv"
# with open(outfile, 'a') as f_handle:
#     np.savetxt(outfile, Ref, delimiter=",")
plt.subplot(211)
n = {0: n_lo, 1: n_hi}  # jakie wycinki dochodza w danej sekundzie czasu
plt.plot(list(map(lambda x: n[x], material)))
plt.xlabel('Time (s)')
plt.ylabel('Refractive index')

plt.subplot(212)
plt.plot(Ref)
plt.ylim([0, 1])
plt.xlabel('Time (s)')
plt.ylabel('Reflectivity')
plt.title('lambda = ' + str(wl) + ' nm')
plt.tight_layout()

# %% rysowanie
root = '/mnt/Dane/0_Microcavities/wneki_symulacje_TMM/refl_vs_growth_t'
params_fig_tgt = root + '/sym_refl_vs_growth_time_' + str(wl)  + '.png'
plt.savefig(params_fig_tgt, dpi=300)
plt.close()
