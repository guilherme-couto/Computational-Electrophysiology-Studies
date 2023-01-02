import matplotlib.pyplot as plt

u = []

with open('minmodel.txt') as f:
    for line in f:
        u.append(float(line))
f.close()

plt.plot(u)
plt.title('Minimal Model for M cells')
plt.ylabel('u (mV)')
plt.xlabel('t (ms)')
plt.savefig('minmodel.png')