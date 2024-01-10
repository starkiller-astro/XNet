'''
This script is used to calculate the molar fraction of C12 and O16 in the
C12 + O16 reaction. It is used to generate CO ratios for the autoencoder.
'''
for i in range(11):
    c12_mf = i * 0.1
    o16_mf = 1 - c12_mf
    c12_Mf = c12_mf / 12
    o16_Mf = o16_mf / 16
    print(f'Run{i + 1}: C12 Molar Fraction = {c12_Mf:.6E}, O16 Molar Fraction = {o16_Mf:.6E}')
    #print("{:.6E}".format(value))
