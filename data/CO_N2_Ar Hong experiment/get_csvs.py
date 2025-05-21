import pandas as pd
import os

# Имя исходного файла
xlsx_file = "Supplementary material S2 raw data.xlsx"

# Загружаем все листы, кроме первого
xls = pd.ExcelFile(xlsx_file)
sheet_names = xls.sheet_names[1:]  # пропускаем первый лист

# Новый заголовок для 10-й строки
new_header = ['Time ms', 'Tvib mean K', 'Trot mean K', 
              'Time ms', 'Ttr cal isentropic K', 'Time ms', 'Pressure bar']

# Обрабатываем каждый лист
for i, sheet_name in enumerate(sheet_names, start=1):
    df = pd.read_excel(xlsx_file, sheet_name=sheet_name, header=None)
    df.replace('--', float('nan'), inplace=True)
    df.dropna(how='all', inplace=True)

    # Создаем директорию для листа
    dir_name = f"Mixture{i}"
    os.makedirs(dir_name, exist_ok=True)

    # Разделяем эксперименты по пустым столбцам
    current_experiment = 1
    start_col = 0
    for col in range(df.shape[1] + 1):
        if col == df.shape[1] or df.iloc[:, col].isnull().all():
            if start_col < col:
                experiment_df = df.iloc[:, start_col:col]
                experiment_df.dropna(how='all', inplace=True)
                experiment_df.dropna(axis=1, how='all', inplace=True)
                if not experiment_df.empty:
                    # Заменить 11-ю строку на новый заголовок
                    experiment_df.iloc[10] = new_header[:experiment_df.shape[1]]
                    
                    # Сохраняем CSV
                    csv_path = os.path.join(dir_name, f"Experiment{current_experiment}.csv")
                    experiment_df.to_csv(csv_path, index=False, header=False)
                    current_experiment += 1
            start_col = col + 1

