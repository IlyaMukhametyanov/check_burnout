import mne
import numpy as np


#legacy
def get_freq(info):
    raw = info
    power, freq = raw.compute_psd().get_data(return_freqs=True)
    return freq



# функция обработки edf > mne.eeg
def input_eeg(filename, type):
    ten_twenty_montage = mne.channels.make_standard_montage('standard_1020')

    if(type =="edf") :
        data = mne.io.read_raw_edf(filename)
        raw_data = data.get_data()
        channels = data.ch_names
        print(channels)
        return data

    elif(type =="bdf"):
        data = mne.io.read_raw_bdf(filename)
        return data


# функция mna.eeg -> ассиметрия альфа ритма
# добавить больше лобных электродов, оформить лучше код
def get_alpha_assimetry(raw, cols):
    freq = get_freq(raw)
    band_names = np.array(['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'Slow', 'Low_beta'])
    filter_list = [[1, 3], [4, 7], [8, 12], [13, 30], [30, 43], [4, 13], [13, 17]]
    bands = []
    for filt in filter_list:
        pt = np.argwhere((freq >= filt[0]) & (freq <= filt[1])).reshape(-1)
        bands.append(pt)
    alpha_frontal = 0
    for index, band in enumerate(bands):
        if (band_names[index] == 'Alpha'):
            power, freq = raw.compute_psd().get_data(return_freqs=True)
            data = power[::, band].mean(axis=1).reshape(1, -1)

            if 'EEG Fp1' in cols and 'EEG Fp2' in cols:
                a_fp1 = data[:, raw.ch_names.index('EEG Fp1')]
                a_fp2 = data[:, raw.ch_names.index('EEG Fp2')]
                beta_fp12 = ((a_fp2 - a_fp1) / (a_fp2 + a_fp1))[0]
                alpha_frontal += beta_fp12

            if 'EEG F7' in cols and 'EEG F8' in cols:
                a_fp1 = data[:, raw.ch_names.index('EEG F7')]
                a_fp2 = data[:, raw.ch_names.index('EEG F8')]
                beta_fp12 = ((a_fp2 - a_fp1) / (a_fp2 + a_fp1))[0]
                alpha_frontal += beta_fp12

            if 'EEG F3' in cols and 'EEG F4' in cols:
                a_fp1 = data[:, raw.ch_names.index('EEG F3')]
                a_fp2 = data[:, raw.ch_names.index('EEG F4')]
                beta_fp12 = ((a_fp2 - a_fp1) / (a_fp2 + a_fp1))[0]
                alpha_frontal += beta_fp12




    return alpha_frontal

def get_beta_assimetry(raw, cols):
    freq = get_freq(raw)
    band_names = np.array(['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'Slow', 'Low_beta'])
    filter_list = [[1, 3], [4, 7], [8, 12], [13, 30], [30, 43], [4, 13], [13, 17]]
    bands = []
    for filt in filter_list:
        pt = np.argwhere((freq >= filt[0]) & (freq <= filt[1])).reshape(-1)
        bands.append(pt)
    beta_assimetry = 0
    count = 0
    hight = 0
    lower = 100
    for index, band in enumerate(bands):
        if (band_names[index] == 'Beta'):
            power, freq = raw.compute_psd().get_data(return_freqs=True)
            data = power[::, band].mean(axis=1).reshape(1, -1)

            if 'EEG Fp1' in cols and 'EEG Fp2' in cols:
                a_fp1 = data[:, raw.ch_names.index('EEG Fp1')]
                a_fp2 = data[:, raw.ch_names.index('EEG Fp2')]
                beta_fp12 = ((a_fp2 - a_fp1) / (a_fp2 + a_fp1))[0]
                beta_assimetry += beta_fp12
                count+=1
                if beta_fp12 > hight:
                    hight = beta_fp12
                if abs(beta_fp12) < lower:
                    lower = abs(beta_fp12)

            if 'EEG F7' in cols and 'EEG F8' in cols:
                a_f7 = data[:, raw.ch_names.index('EEG F7')]
                a_f8 = data[:, raw.ch_names.index('EEG F8')]
                beta_f78 = ((a_f7 - a_f8) / (a_f7 + a_f8))[0]
                beta_assimetry += beta_f78
                count+=1
                if beta_f78 > hight:
                    hight = beta_f78
                if abs(beta_f78) < lower:
                    lower = abs(beta_f78)

            if 'EEG F3' in cols and 'EEG F4' in cols:
                a_f3 = data[:, raw.ch_names.index('EEG F3')]
                a_f4 = data[:, raw.ch_names.index('EEG F4')]
                beta_f34 = ((a_f3 - a_f4) / (a_f3 + a_f4))[0]
                beta_assimetry += beta_f34
                count+=1
                if beta_f34 > hight:
                    hight = beta_f34
                if abs(beta_f34) < lower:
                    lower = abs(beta_f34)

            if 'EEG T3' in cols and 'EEG T4' in cols:
                a_t3 = data[:, raw.ch_names.index('EEG T3')]
                a_t4 = data[:, raw.ch_names.index('EEG T4')]
                beta_t34 = ((a_t3 - a_t4) / (a_t3 + a_t4))[0]
                beta_assimetry += beta_t34
                count+=1
                if beta_t34 > hight:
                    hight = beta_t34
                if abs(beta_t34) < lower:
                    lower = abs(beta_t34)

            if 'EEG C3' in cols and 'EEG C4' in cols:
                a_c3 = data[:, raw.ch_names.index('EEG C3')]
                a_c4 = data[:, raw.ch_names.index('EEG C4')]
                beta_c34 = ((a_c3 - a_c4) / (a_c3 + a_c4))[0]
                beta_assimetry += beta_c34
                count+=1
                if beta_c34 > hight:
                    hight = beta_c34
                if abs(beta_c34) < lower:
                    lower = abs(beta_c34)

            if 'EEG P3' in cols and 'EEG P4' in cols:
                a_p3 = data[:, raw.ch_names.index('EEG P3')]
                a_p4 = data[:, raw.ch_names.index('EEG P4')]
                beta_p34 = ((a_p3 - a_p4) / (a_p3 + a_p4))[0]
                beta_assimetry += beta_p34
                count+=1
                if beta_p34 > hight:
                    hight = beta_p34
                if abs(beta_p34) < lower:
                    lower = abs(beta_p34)

            if 'EEG T5' in cols and 'EEG T6' in cols:
                a_t5 = data[:, raw.ch_names.index('EEG T5')]
                a_t6 = data[:, raw.ch_names.index('EEG T6')]
                beta_t56 = ((a_t5 - a_t6) / (a_t5 + a_t6))[0]
                beta_assimetry += beta_t56
                count+=1
                if beta_t56 > hight:
                    hight = beta_t56
                if abs(beta_t56) < lower:
                    lower = abs(beta_t56)

            if 'EEG O1' in cols and 'EEG O2' in cols:
                a_o1 = data[:, raw.ch_names.index('EEG O1')]
                a_o2 = data[:, raw.ch_names.index('EEG O2')]
                beta_012 = ((a_o1 - a_o2) / (a_o1 + a_o2))[0]
                beta_assimetry += beta_012
                count+=1
                if beta_012 > hight:
                    hight = beta_012
                if abs(beta_012) < lower:
                    lower = abs(beta_012)

    return (beta_assimetry/ count), hight, lower
#функции оценки состояния человека

#функции вычисления альфа,бета
def get_alpha_beta_average(raw, col):
    freq = get_freq(raw)
    band_names = np.array(['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'Slow', 'Low_beta'])
    filter_list = [[1, 3], [4, 7], [8, 12], [13, 30], [30, 43], [4, 13], [13, 17]]
    bands = []
    for filt in filter_list:
        pt = np.argwhere((freq >= filt[0]) & (freq <= filt[1])).reshape(-1)
        bands.append(pt)

    one_alpha = []
    one_beta = []

    for index, band in enumerate(bands):
        if (band_names[index] == 'Alpha'):
            power, freq = raw.compute_psd().get_data(return_freqs=True)
            data = power[::, band].mean(axis=1).reshape(1, -1)

            for x in col:
                item_1 = data[:, raw.ch_names.index(x)]
                one_alpha.append(item_1)
        if (band_names[index] == 'Beta'):
            power, freq = raw.compute_psd().get_data(return_freqs=True)
            data = power[::, band].mean(axis=1).reshape(1, -1)

            for x in col:
                item_2 = data[:, raw.ch_names.index(x)]
                one_beta.append(item_2)


    averge_dif_alpha_subject = sum(one_alpha) / len(one_alpha)
    averge_dif_bets_subject = sum(one_beta) / len(one_beta)
    return averge_dif_alpha_subject, averge_dif_bets_subject

#функции выявления паттернов выгорания

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    filename = "real_data\\files\\S001\\S001R04.edf"

    raw = input_eeg(filename, "edf")
    channels = raw.ch_names
    print(channels)

    # ['Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'T3', 'T4', 'Fz', 'Cz', 'Pz', 'Aux']
    # ['EEG Fp1', 'EEG Fp2', 'EEG F3', 'EEG F4', 'EEG F7', 'EEG F8', 'EEG T3', 'EEG T4', 'EEG C3', 'EEG C4', 'EEG T5',
    #  'EEG T6', 'EEG P3', 'EEG P4', 'EEG O1', 'EEG O2', 'EEG Fz', 'EEG Cz', 'EEG Pz', 'EEG A2-A1', 'ECG ECG']

    # file - EEG Motor Movement - ['Fc5.', 'Fc3.', 'Fc1.', 'Fcz.', 'Fc2.', 'Fc4.', 'Fc6.', 'C5..', 'C3..', 'C1..', 'Cz..', 'C2..', 'C4..', 'C6..',
    #  'Cp5.', 'Cp3.', 'Cp1.', 'Cpz.', 'Cp2.', 'Cp4.', 'Cp6.', 'Fp1.', 'Fpz.', 'Fp2.', 'Af7.', 'Af3.', 'Afz.', 'Af4.',
    #  'Af8.', 'F7..', 'F5..', 'F3..', 'F1..', 'Fz..', 'F2..', 'F4..', 'F6..', 'F8..', 'Ft7.', 'Ft8.', 'T7..', 'T8..',
    #  'T9..', 'T10.', 'Tp7.', 'Tp8.', 'P7..', 'P5..', 'P3..', 'P1..', 'Pz..', 'P2..', 'P4..', 'P6..', 'P8..', 'Po7.',
    #  'Po3.', 'Poz.', 'Po4.', 'Po8.', 'O1..', 'Oz..', 'O2..', 'Iz..']

    # alpha_frontal = get_alpha_assimetry(raw, channels)
    # print("Ассиметрия альфа ритма в лобной доли ", alpha_frontal)
    #
    # averge_dif_alpha_subject, averge_dif_bets_subject = get_alpha_beta_average(raw,channels)
    # print("Средние альфа и бета ритмы ",averge_dif_alpha_subject, averge_dif_bets_subject)
    #
    # beta_assimenry, hight, lower = get_beta_assimetry(raw,channels)
    # print("Межполушарная ассиметрия бета ритма, средняя, наибольша, наименьшая")
    # print(beta_assimenry, hight, lower )

# See PyCharm help at https://www.jetbrains.com/help/pycharm/