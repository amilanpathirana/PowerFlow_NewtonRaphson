import numpy as np
from createbus import Bus
from createline import Line

# Base Power
Sb = 100.0

# NR variables
tolerance = 0.00001
counter = 0


# read data from file and split them properly

readData = open("inputdata.txt", "r").read()

dataSection = readData.split("----\n")
print('dataSection[0]')
print(dataSection[0])
print('dataSection[1]')
print(dataSection[1])
# Create bus objects
busdata = dataSection[0].split('\n')
print('busdata[0]')
print(busdata[0])
print('busdata[1]')
print(busdata[1])
print('busdata[2]')
print(busdata[2])

buses = dict()
for item in busdata:
    if item.strip():
        buses[str(int(item[0:4]))] = Bus(item, Sb)

# Create line objects
lines = dict()
linedata = dataSection[1].split("\n")

for item in linedata:
    if item.strip():
        lines[item[0:4].strip() + "-" + item[4:12].strip()] = Line(item, Sb)
print(lines)
print('///////////////////////////////////////////')
print(buses)


# Create file
f = open("results.txt", "w+")

# Nodal Admittance Matrix
Ybus = np.zeros((len(buses), len(buses)), dtype=complex)

# Shunt Elements Vector
Bshunt = np.zeros(len(buses), dtype=complex)

for key in lines.keys():
    Ybus[lines[key].origin - 1][lines[key].destiny - 1] = -1/(lines[key].R + 1j*lines[key].X)
    Bshunt[lines[key].origin - 1] += 1j*lines[key].B/2
    Bshunt[lines[key].destiny - 1] += 1j*lines[key].B/2

Ybus += Ybus.T

np.fill_diagonal(Ybus, Bshunt - np.sum(Ybus, axis=1))

# Flat Start
for key in buses.keys():
    if buses[key].bustype == 'PQ':
        buses[key].V = 1.
        buses[key].theta = 0.
    elif buses[key].bustype == 'PV':
        buses[key].theta = 0.

# Initialize Active and Reactive Power
P = np.zeros(len(buses))
Q = np.zeros(len(buses))

# Initialize Mismatches
misP = np.zeros(len(buses))
misQ = np.zeros(len(buses))

# Power Calculation
for bus in range(len(buses)):
    for otherbus in range(len(buses)):

        # Calculate angle difference
        theta_km = buses[str(bus + 1)].theta - buses[str(otherbus + 1)].theta

        # Calculate active and reactive power reaching bus
        P[bus] += buses[str(bus+1)].V*buses[str(otherbus+1)].V*(np.real(Ybus[bus, otherbus])*np.cos(theta_km) + np.imag(Ybus[bus, otherbus]*np.sin(theta_km)))
        Q[bus] += buses[str(bus+1)].V*buses[str(otherbus+1)].V*(np.real(Ybus[bus, otherbus])*np.sin(theta_km) - np.imag(Ybus[bus, otherbus]*np.cos(theta_km)))

    # Save calculated values of power
    buses[str(bus+1)].save_power(P[bus], Q[bus])

    # Calculate mismatches
    if buses[str(bus+1)].bustype == 'PQ':
        misP[bus] = buses[str(bus+1)].P - P[bus]
        misQ[bus] = buses[str(bus+1)].Q - Q[bus]
    elif buses[str(bus+1)].bustype == 'PV':
        misP[bus] = buses[str(bus+1)].P - P[bus]

# Mismatch vector
mis = np.vstack((np.array([misP]).T, np.array([misQ]).T))

# Print Status
print("Iteration #%d" % counter)
print("Error: %f" % max(abs(mis))[0])
for key in buses.keys():
    print("V%s = %4f < %2f" % (key, buses[key].V, buses[key].theta))
print(30*"-")

# Write in results file
f.write("Iteration #%d" % counter)
f.write("\nError: %f" % max(abs(mis))[0])
for key in buses.keys():
    f.write("\nV%s = %4f < %2f" % (key, buses[key].V, buses[key].theta))
f.write("\n" + 30*"-" + "\n")

while max(abs(mis)) > tolerance and counter < 100:

    # Create Jacobian submatrices
    H = np.zeros((len(buses), len(buses)))
    N = np.zeros((len(buses), len(buses)))
    M = np.zeros((len(buses), len(buses)))
    L = np.zeros((len(buses), len(buses)))

    # Calculate Jacobian submatrices
    for bus in range(len(buses)):
        for otherbus in range(len(buses)):

            # Calculate angle difference
            theta_km = buses[str(bus + 1)].theta - buses[str(otherbus + 1)].theta

            # Calculate values outside main diagonal
            if bus is not otherbus:
                H[bus, otherbus] = buses[str(bus+1)].V*buses[str(otherbus+1)].V*(np.real(Ybus[bus, otherbus])*np.sin(theta_km) - np.imag(Ybus[bus, otherbus]*np.cos(theta_km)))
                N[bus, otherbus] = buses[str(bus+1)].V*(np.real(Ybus[bus, otherbus])*np.cos(theta_km) + np.imag(Ybus[bus, otherbus]*np.sin(theta_km)))
                M[bus, otherbus] = -buses[str(bus+1)].V*buses[str(otherbus+1)].V*(np.real(Ybus[bus, otherbus])*np.cos(theta_km) + np.imag(Ybus[bus, otherbus]*np.sin(theta_km)))
                L[bus, otherbus] = buses[str(bus+1)].V*(np.real(Ybus[bus, otherbus])*np.sin(theta_km) - np.imag(Ybus[bus, otherbus]*np.cos(theta_km)))
            # Calculate values in main diagonal
            else:
                # Don't apply changes of angle to Vθ buses
                if buses[str(bus+1)].bustype == 'Vθ':
                    H[bus, otherbus] = 10**99
                else:
                    H[bus, otherbus] = -Q[bus] - (buses[str(bus+1)].V**2)*np.imag(Ybus[bus, bus])

                N[bus, otherbus] = (P[bus] + (buses[str(bus+1)].V**2)*np.real(Ybus[bus, bus]))/buses[str(bus+1)].V
                M[bus, otherbus] = P[bus] - (buses[str(bus+1)].V**2)*np.real(Ybus[bus, bus])

                # Don't apply changes of voltage to Vθ and PV buses
                if buses[str(bus+1)].bustype == 'PQ':
                    L[bus, otherbus] = (Q[bus] - (buses[str(bus+1)].V**2)*np.imag(Ybus[bus, bus]))/buses[str(bus+1)].V
                else:
                    L[bus, otherbus] = 10**99

    J = np.vstack((np.hstack((H, N)), np.hstack((M, L))))

    correction = np.linalg.solve(J, mis)

    for index in range(int(len(correction)/2)):
        buses[str(index + 1)].refresh(correction[index][0], correction[index + len(buses)][0])

    # Initialize Active and Reactive Power
    P = np.zeros(len(buses))
    Q = np.zeros(len(buses))

    # Initialize Mismatches
    misP = np.zeros(len(buses))
    misQ = np.zeros(len(buses))

    # Power Calculation
    for bus in range(len(buses)):
        for otherbus in range(len(buses)):

            # Calculate angle difference
            theta_km = buses[str(bus + 1)].theta - buses[str(otherbus + 1)].theta

            # Calculate active and reactive power reaching bus
            P[bus] += buses[str(bus+1)].V*buses[str(otherbus+1)].V*(np.real(Ybus[bus, otherbus])*np.cos(theta_km) + np.imag(Ybus[bus, otherbus]*np.sin(theta_km)))
            Q[bus] += buses[str(bus+1)].V*buses[str(otherbus+1)].V*(np.real(Ybus[bus, otherbus])*np.sin(theta_km) - np.imag(Ybus[bus, otherbus]*np.cos(theta_km)))

        # Save calculated values of power
        buses[str(bus + 1)].save_power(P[bus], Q[bus])

        # Calculate mismatches
        if buses[str(bus+1)].bustype == 'PQ':
            misP[bus] = buses[str(bus+1)].P - P[bus]
            misQ[bus] = buses[str(bus+1)].Q - Q[bus]
        elif buses[str(bus+1)].bustype == 'PV':
            misP[bus] = buses[str(bus+1)].P - P[bus]

    # Mismatch vector
    mis = np.vstack((np.array([misP]).T, np.array([misQ]).T))

    # Refresh counter
    counter += 1

    # Print Status
    print("\nIteration #%d" % counter)
    print("Error: %f" % max(abs(mis))[0])
    for key in buses.keys():
        print("V%s = %4f < %2f" % (key, buses[key].V, buses[key].theta))
    print(30*"-")

    # Write in results file
    f.write("\nIteration #%d" % counter)
    f.write("\nError: %f" % max(abs(mis))[0])
    for key in buses.keys():
        f.write("\nV%s = %4f < %2f" % (key, buses[key].V, buses[key].theta))
    f.write("\n" + 30*"-" + "\n")

# Power flow calculation saving it on line object
for key in lines.keys():

    theta_km = buses[str(lines[key].origin)].theta - buses[str(lines[key].destiny)].theta

    Vk = buses[str(lines[key].origin)].V
    Vm = buses[str(lines[key].destiny)].V

    Y = 1/(lines[key].R + 1j*lines[key].X)

    pkm = (Vk**2)*np.real(Y) - Vk*Vm*(np.real(Y)*np.cos(theta_km) + np.imag(Y)*np.sin(theta_km))
    qkm = -(Vk**2)*(np.imag(Y) + lines[key].B/2) - Vk*Vm*(np.real(Y)*np.sin(theta_km) - np.imag(Y)*np.cos(theta_km))

    lines[key].save_flow(pkm, qkm, lines[key].origin)

    pmk = (Vm**2)*np.real(Y) - Vm*Vk*(np.real(Y)*np.cos(-theta_km) + np.imag(Y)*np.sin(-theta_km))
    qmk = -(Vm**2)*(np.imag(Y) + lines[key].B/2) - Vk*Vm*(np.real(Y)*np.sin(-theta_km) - np.imag(Y)*np.cos(-theta_km))

    lines[key].save_flow(pmk, qmk, lines[key].destiny)

# Print power flow results
print("\nPower Flow\n")
print("Line", 6*"\t", "Skm", 11*"\t", "Smk", 11*"\t", "Perdas")
print(150*"-")

for key in lines.keys():
    print(key, "\t", lines[key].S_od, "\t", lines[key].S_do, "\t", lines[key].S_od - lines[key].S_do)

# Write power flow results on file
f.write("\n\nPower Flow\n")
f.write("\n" + "Line" + 6*"\t" + "Skm" + 11*"\t" + "Smk" + 11*"\t" + "Perdas")
f.write("\n" + 150*"-")
for key in lines.keys():
    f.write("\n" + key + "\t" + str(lines[key].S_od) + "\t" + str(lines[key].S_do) + "\t" + str(lines[key].S_od - lines[key].S_do))
