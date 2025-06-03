# paperSO-SSM
Code the numerical experiments of the paper "Self-organizing state-space models with artifical dymanics" (available at https://arxiv.org/abs/2409.08928)

Warning: the dataset for the expriments in Section 4.1.2 is not open access. In this repository this dataset has been replaced by observations for the city of Portland (US),  available on Kaggle (at https://www.kaggle.com/datasets/selfishgene/historical-hourly-weather-data?select=temperature.csv_). The results are not as good as in the paper, which can be due (i) to a difference in the quality of the data, (ii) to the fact that q=4 basis functions is not enough for this dataset or (iii) to the fact that our splitting of each year in two six-month period is not appropriate for the temperature in Portland.
