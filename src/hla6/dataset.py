import os
import sys
import pandas as pd
import gzip
import yaml
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

class CustomDataset(Dataset):
    def __init__(self, features_file='features.txt', labels_file='labels.txt', label_encoding='label_encoding.txt', out_file='hla'):
        self.out_file = out_file
        df = pd.read_table(features_file, header=0, sep='\t')
        mat = df.iloc[:, 1:].values
        mat = mat.reshape(mat.shape[0], mat.shape[1]//2, 2).astype('float32')
        self.X = torch.tensor(mat).permute(0, 2, 1)

        df = pd.read_table(labels_file, header=0, sep='\t')
        mat = df.iloc[:, 1:].values
        mat = mat.reshape(mat.shape[0], mat.shape[1]//2, 2).astype('int')
        self.y = torch.tensor(mat).permute(0, 2, 1)
        print(self.X.shape)
        print(self.y.shape)

    def __len__(self):
        return(len(self.y))

    def __getitem__(self, idx):
        X = self.X[idx]
        y = self.y[idx]
        return(X, y)

    def split_save_dataset(self, ratio=[0.8, 0.1, 0.1], batch_size=32, shuffle=True, num_workers=0):
        self.train_file = self.out_file + '_dataset_train.pt'
        self.val_file = self.out_file + '_dataset_val.pt'
        self.test_file = self.out_file + '_dataset_test.pt'

        self.ds_train, self.ds_val, self.ds_test = torch.utils.data.random_split(self, ratio)

        self.dl_train = DataLoader(self.ds_train, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
        self.dl_val = DataLoader(self.ds_val, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
        self.dl_test = DataLoader(self.ds_test, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
        torch.save(self.dl_train, self.train_file)
        torch.save(self.dl_val, self.val_file)
        torch.save(self.dl_test, self.test_file)

if __name__ == '__main__':
    ds = CustomDataset()
    ds.split_save_dataset()

