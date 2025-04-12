# Antoni Rumowski 193087 Marek Marcinko 197870

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import control as ct

print("Imports complete.")
print("NumPy version:", np.__version__)
print("Matplotlib version:", matplotlib.__version__)
print("Control version:", ct.__version__)
input("Press Enter to continue...")

is_stable = True

# input coefficients

tran_numerator = float(input("b0: "))
a3  = float(input("a3: "))
a2 = float(input("a2: "))
a1 = float(input("a1: "))
a0 = float(input("a0: "))
print("s**4 +", a3, "* s**3+", a2, "* s**2+", a1, "* s+", a0)
tran_denominator = [1, a3, a2, a1, a0]

def bode_plot(G):
    mag, phase, omega = ct.bode(G, plot=False)
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.semilogx(omega, 20 * np.log10(mag))  # log scale for frequency
    plt.title('Bode Plot')
    plt.ylabel('Magnitude (dB)')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.semilogx(omega, np.degrees(phase))
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('Phase (deg)')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig("bode_plot.png")
    plt.close()
    print("Bode plot complete.")

def step_response(G):
    t = np.linspace(0, 10, 1000)  # time
    t_out, y_out = ct.step_response(G, T=t)

    plt.figure()
    plt.plot(t_out, y_out, label="y(t) – step response")
    plt.plot(t_out, np.ones_like(t_out), '--', label="u(t) – unit step")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid()
    plt.title("Step response of the system")
    plt.savefig("step_response.png")
    plt.close()
    print("Step response complete.")

def sine_response(G):
    w = 1.0  #sine frequency
    t = np.linspace(0, 20, 1000)
    u = np.sin(w * t)

    t_out, y_out = ct.forced_response(G, T=t, U=u)

    plt.figure()
    plt.plot(t_out, u, '--', label="u(t) – sine input")
    plt.plot(t_out, y_out, label="y(t) – system response")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid()
    plt.title("System response to sine input")
    plt.savefig("sine_response.png")
    plt.close()
    print("Sine response complete.")

def check_stability(tab):
    rows2 = [0,1]
    cols2 = [0,1]
    matrix2 = tab[np.ix_(rows2, cols2)]

    test2 = np.linalg.det(matrix2)

    rows3 = [0,1,2]
    cols3 = [0,1,2]
    matrix3 = tab[np.ix_(rows3, cols3)]

    test3 = np.linalg.det(matrix3)
    test4 = np.linalg.det(tab)

    print("Matrix 2x2:", matrix2)
    print("Matrix 3x3:", matrix3)
    if test2 <= 0 or test3 <= 0 or test4 <= 0:
        return False
    else:
        return True
    
#check for negative coefficients in the characteristic polynomial

for nums in tran_denominator:
    if nums < 0:
        is_stable = False
        print("The system is unstable, cannot simulate the output.")
        break

tab_hurwitz = np.array([[a3,a1,0,0],
                     [1,a2,a0,0],
                     [0,a3,a1,0],
                     [0,1,a2,a0]])
print(tab_hurwitz)
is_stable = check_stability(tab_hurwitz)
G = ct.TransferFunction([tran_numerator], tran_denominator)
print("Transfer function:", G)

#the result of the Hurwitz test
if is_stable:
    print("The system is stable.")
else:   
    print("The system is unstable.")

bode_plot(G)
input("Press Enter to continue...") 
step_response(G)
input("Press Enter to continue...")
sine_response(G)
print("All tasks complete.")
input("Press Enter to exit...")
SystemExit(0)


