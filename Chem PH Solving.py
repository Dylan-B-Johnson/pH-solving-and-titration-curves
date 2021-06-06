#    Copyright [2021] [Dylan Johnson]

#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at

#        http://www.apache.org/licenses/LICENSE-2.0

#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.


import matplotlib.pyplot as plt
import math

# At the equivlance point (when the analyte has been completely neutralized), the pH of the solution is determined
# by the ionization of the salt in water
# E.G. pH of CH3COOH + NaOH -> CH3COO- Na+ + HOH @ equivlance is determined by:
# CH3COO- + HOH <--> CH3COOH + OH-  (here Kb=[CH3COOH][OH-]/[CH3COO-] and Kb is calculated from the Ka of CH3COOH)
# Or HCl + NaOH -> NaCl + HOH @ equivlance is pH=7 (NaCl is strong base / strong acid; therefore neutral)

log = lambda x: math.log(x, 10)
switch = lambda k: 10 ** (-(14 - p(k)))
p = lambda x: -log(x)
p_quad = lambda a, b, c: ((-b + (b * b - 4 * a * c) ** 0.5) / (2 * a))
n_quad = lambda a, b, c: ((-b - (b * b - 4 * a * c) ** 0.5) / (2 * a))


# --------------react----------------------
# Solves Stoich neutralizations in the form: (Can be used for any reaction with 2 reactants and 1-2 products)
# Analyte + Titrant --> Salt + Water
# Or 
# Analyte + Titrant --> Salt
# ratio = [1,1,1] or [1,1,1,1] 
# NOTE: ALWAYS ASSUMES
# C1 -> [Analyte]
# C2 -> [Titrant]
# V1 -> Volume of Analyte
# V2 -> Volume of Titrant
# Unit -> "mL" or "L"
# Returns {'analyte': mols of Analyte remaining, 'titrant': mols of Titrant remaining, 'salt': mols of Salt remaining,
# 'V': volume of solution in L, 'titrant_mol_needed': moles of Titrant needed to reach equivlance pt,
# 'titrant_mol_needed': volume (in L) of Titrant needed to reach equivlance pt}
def react(ratio=[1, 1, 1, 1], C1=5, C2=5, V1=50, V2=25, unit='mL'):
    if unit == 'mL':
        mol1 = (V1 / 1000) * C1
        mol2 = (V2 / 1000) * C2
        V = V1 / 1000 + V2 / 1000
    elif unit == 'L':
        mol1 = V1 * C1
        mol2 = V2 * C2
        V = V1 + V2
    else:
        print('Unit Error')
    titrant_mol_needed = mol1 / ratio[0] * ratio[1]
    titrant_vol_needed = titrant_mol_needed / C2
    if len(ratio) == 3 or len(ratio) == 4:
        if len(ratio) == 3:
            water = 0.0
        if (mol1 / ratio[0] <= mol2 / ratio[1]):
            analyte = 0.0
            titrant = mol2 - mol1 / ratio[0] * ratio[1]
            salt = mol1 / ratio[0] * ratio[2]
            if len(ratio) == 4:
                water = mol1 / ratio[0] * ratio[3]
            return {'analyte': analyte, 'titrant': titrant, 'salt': salt, 'V': V, 'water': water,
                    'titrant_mol_needed': titrant_mol_needed, 'titrant_vol_needed': titrant_vol_needed}
        if (mol1 / ratio[0] > mol2 / ratio[1]):
            analyte = mol1 - mol2 / ratio[1] * ratio[0]
            titrant = 0.0
            salt = mol2 / ratio[1] * ratio[2]
            if len(ratio) == 4:
                water = mol2 / ratio[1] * ratio[3]
            return {'analyte': analyte, 'titrant': titrant, 'salt': salt, 'V': V, 'water': water,
                    'titrant_mol_needed': titrant_mol_needed, 'titrant_vol_needed': titrant_vol_needed}
    else:
        print('Ratio Length Error')


# acid_or_base 'acid' for buffers where the analyte is an acid, and 'base' for buffers where the analyte is a base
# ka if 'acid', kb if 'base'
def buffer_pH(acid_or_base, k, analyte_mol, salt_mol):
    if acid_or_base == 'acid':
        return (p(k) + log(salt_mol / analyte_mol))
    if acid_or_base == 'base':
        return 14 - ((p(k) + log(salt_mol / analyte_mol)))


# acid_or_base = 'acid' or 'base'
# TAKES CONCENTRATION, NOT MOLES
def strong_pH(acid_or_base, conc):
    if acid_or_base == 'acid':
        return (p(conc))
    if acid_or_base == 'base':
        return (14 - p(conc))


# acid_or_base = 'acid' or 'base'
# TAKES CONCENTRATION, NOT MOLES
def weak_pH(acid_or_base, k, conc):
    if abs(p_quad(-1, -k, k * conc)) != abs(n_quad(-1, -k, k * conc)):
        #         print('Error: Quadratic did not yield consistent X value.')
        if acid_or_base == 'acid':
            H = (k * conc) ** 0.5
            return (p(H))
        if acid_or_base == 'base':
            H = switch((k * conc) ** 0.5)
            return (p(H))
    else:
        if acid_or_base == 'acid':
            H = abs(p_quad(-1, -k, k * conc))
            return (p(H))
        if acid_or_base == 'base':
            H = switch(abs(p_quad(-1, -k, k * conc)))
            return (p(H))


# TAKES CONCENTRATION, NOT MOLES
# acid_or_base is for analyte 
def equivlance_pH(acid_or_base, k, salt_conc):
    # weak_pH but convert K to opposite
    if acid_or_base == 'acid':
        return weak_pH('base', switch(k), salt_conc)
    if acid_or_base == 'base':
        return weak_pH('acid', switch(k), salt_conc)


# rxn_mol -> output of react(), or:
#      {'analyte': mols of Analyte remaining, 'titrant': mols of Titrant remaining, 'salt': mols of Salt remaining,
# 'V': volume of solution in L, 'titrant_mol_needed': moles of Titrant needed to reach equivlance pt,
# 'titrant_mol_needed': volume (in L) of Titrant needed to reach equivlance pt}
# Ka/Kb -> constant for analyte (0 for the one that does not apply) (unless strong analyte)
# acid_or_base is for analyte ('acid' or 'base')
# if strong analyte - weak titrant, k is the k of the titrant
# k2 is only for weak-weak titrations, where it is the k of the titrant 
def get_pH(rxn_mol, k=1.7e-5, acid_or_base='acid', strong_titrant=True, strong_analyte=False, k2=1.8e-5):
    # weak analyte, strong titrant
    if not (strong_analyte) and strong_titrant:
        # initial pH
        if rxn_mol['analyte'] > 0 and rxn_mol['titrant'] == 0 and rxn_mol['salt'] == 0:
            return weak_pH(acid_or_base, k, rxn_mol['analyte'] / rxn_mol['V'])
        # before equivlance pt (buffer solution)
        if rxn_mol['analyte'] > 0 and rxn_mol['titrant'] == 0:
            if acid_or_base == 'acid':
                return buffer_pH(acid_or_base, k, rxn_mol['analyte'], rxn_mol['salt'])
            if acid_or_base == 'base':
                return buffer_pH(acid_or_base, k, rxn_mol['analyte'], rxn_mol['salt'])
        # @ equivlance pt
        if rxn_mol['analyte'] == 0 and rxn_mol['titrant'] == 0:
            return equivlance_pH(acid_or_base, k, (rxn_mol['salt'] / rxn_mol['V']))
        # past equivlance pt (strong acid/base for this combo)
        if rxn_mol['titrant'] > 0 and rxn_mol['analyte'] == 0:
            if acid_or_base == 'acid':
                return strong_pH('base', rxn_mol['titrant'] / rxn_mol['V'])
            if acid_or_base == 'base':
                return strong_pH('acid', rxn_mol['titrant'] / rxn_mol['V'])
    elif strong_analyte and not (strong_titrant):
        # initial pH
        if rxn_mol['analyte'] > 0 and rxn_mol['titrant'] == 0 and rxn_mol['salt'] == 0:
            return strong_pH(acid_or_base, rxn_mol['analyte'] / rxn_mol['V'])
        # before equivlance pt (buffer solution)
        if rxn_mol['analyte'] > 0 and rxn_mol['titrant'] == 0:
            return strong_pH(acid_or_base, rxn_mol['analyte'] / rxn_mol['V'])
        # @ equivlance pt
        if rxn_mol['analyte'] == 0 and rxn_mol['titrant'] == 0:
            return 7.0
        # past equivlance pt (strong acid/base for this combo)
        if rxn_mol['titrant'] > 0 and rxn_mol['analyte'] == 0:
            if acid_or_base == 'acid':
                return weak_pH('base', k, rxn_mol['titrant'] / rxn_mol['V'])
            if acid_or_base == 'base':
                return weak_pH('acid', k, rxn_mol['titrant'] / rxn_mol['V'])
    elif strong_analyte and strong_titrant:
        # initial pH
        if rxn_mol['analyte'] > 0 and rxn_mol['titrant'] == 0 and rxn_mol['salt'] == 0:
            return strong_pH(acid_or_base, rxn_mol['analyte'] / rxn_mol['V'])
        # before equivlance pt (buffer solution)
        if rxn_mol['analyte'] > 0 and rxn_mol['titrant'] == 0:
            return strong_pH(acid_or_base, rxn_mol['analyte'] / rxn_mol['V'])
        # @ equivlance pt
        if rxn_mol['analyte'] == 0 and rxn_mol['titrant'] == 0:
            return 7.0
        # past equivlance pt (strong acid/base for this combo)
        if rxn_mol['titrant'] > 0 and rxn_mol['analyte'] == 0:
            if acid_or_base == 'acid':
                return strong_pH('base', rxn_mol['titrant'] / rxn_mol['V'])
            if acid_or_base == 'base':
                return strong_pH('acid', rxn_mol['titrant'] / rxn_mol['V'])
    elif not (strong_analyte) and not (strong_titrant):
        print('Weak-Weak titrations have not been implimented yet. They are not good experimental design.')


#         # initial pH
#         if rxn_mol['analyte']>0 and rxn_mol['titrant']==0 and rxn_mol['salt']==0:
#             return weak_pH(acid_or_base,k,rxn_mol['analyte']/rxn_mol['V'])
#         # before equivlance pt (buffer solution)
#         if rxn_mol['analyte']>0 and rxn_mol['titrant']==0:
#             if acid_or_base=='acid':
#                 return buffer_pH(acid_or_base,k,rxn_mol['analyte'],rxn_mol['salt'])
#             if acid_or_base=='base':
#                 return buffer_pH(acid_or_base,k,rxn_mol['analyte'],rxn_mol['salt'])
#         # @ equivlance pt
#         if rxn_mol['analyte']==0 and rxn_mol['titrant']==0:
#             return 7.0
#         # past equivlance pt (strong acid/base for this combo)
#         if rxn_mol['titrant']>0 and rxn_mol['analyte']==0:
#             if acid_or_base=='acid':
#                 return strong_pH('base',rxn_mol['titrant']/rxn_mol['V'])
#             if acid_or_base=='base':
#                 return strong_pH('acid',rxn_mol['titrant']/rxn_mol['V'])
# equivlance point of strong-strong if 1:1:1:1 is always pH=7.00
# weak-weak is not a good experimental design

# final_vol is max ammount added 
# if strong analyte - weak titrant, k is the k of the titrant
# k2 is only for weak-weak titrations, where it is the k of the titrant 
def plot_titration(initial_vol=0, final_vol=100, increment=0.1, ratio=[1, 1, 1, 1], C1=5, C2=5, V1=50, V2=25, unit='mL',
                   strong_titrant=True, strong_analyte=False, k=1.7e-5, acid_or_base='acid', k2=1.8e-5):
    added = initial_vol
    x = []
    y = []
    while added <= final_vol:
        y.append(get_pH(react(ratio=ratio, C1=C1, C2=C2, V1=V1, V2=added, unit=unit),
                        k=k, acid_or_base=acid_or_base, strong_titrant=strong_titrant, strong_analyte=strong_analyte))
        x.append(added)
        if y[len(y) - 1] < 0 or y[len(y) - 1] > 14:
            del (y[len(y) - 1])
            del (x[len(x) - 1])
        added += increment
    plt.xlabel('Volume Added ' + '(' + unit + ')')
    plt.ylabel('pH')
    plt.ylim(bottom=0, top=14)
    plt.plot(x, y)
    if acid_or_base == 'acid' and strong_titrant == True and strong_analyte == False:
        plt.title('Weak Acid-Strong Base Titration Curve')
    if acid_or_base == 'base' and strong_titrant == True and strong_analyte == False:
        plt.title('Weak Base-Strong Acid Titration Curve')
    plt.show()
    tit_mol_needed = react(ratio=ratio, C1=C1, C2=C2, V1=V1, V2=V2, unit=unit)['titrant_mol_needed']
    tit_vol_needed = react(ratio=ratio, C1=C1, C2=C2, V1=V1, V2=V2, unit=unit)['titrant_vol_needed'] * 1000
    if not (strong_analyte) and strong_titrant:
        salt_conc_at_eqvl = react(ratio=ratio, C1=C1, C2=C2, V1=V1, V2=tit_vol_needed, unit='mL')['salt'] / \
                            react(ratio=ratio, C1=C1, C2=C2, V1=V1, V2=tit_vol_needed, unit='mL')['V']
        equ_pH = equivlance_pH(acid_or_base, k, salt_conc_at_eqvl)
    elif (strong_analyte and not (strong_titrant)) or (strong_analyte and strong_titrant):
        equ_pH = 7.0
    print('Initial pH:' + str(y[0]))
    print('pH at Equivlance Point: ' + str(equ_pH))
    print('Final pH:', str(y[len(y) - 1]))
    print('Volume of Titrant Needed for Equivlance:', str(tit_vol_needed), 'mL')
    print('Ammount of Titrant Needed for Equivlance:', str(tit_mol_needed), 'mol')
    return x, y


plot_titration()
