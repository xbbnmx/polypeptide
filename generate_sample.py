import pandas as pd
import math
import argparse
from featured_data_generated import cal_pep_des
import time

class GenerateSample():
    def __init__(self, grampa_path, negative_file_path, generate_example_path, mode):
        self.grampa_path = grampa_path
        self.negative_file_path = negative_file_path
        self.generate_example_path = generate_example_path
        self.mode = mode

    def __call__(self, *args, **kwargs):
        data = pd.read_csv(self.grampa_path, encoding="utf8")#读取grampa数据
        data = self.filter_data_with_str("bacterium", "aureus", data)#按照这两列取数
        positive_sample = self.generate_all_peptides(data)#计算数据的sequence，MIC(抑菌浓度),type列
        negtive_sample = self.generate_negative_data(self.negative_file_path)#只取negative文件中sequence列不包含"B|X|Z|O|U"的数据
        all_sample = self.concat_datasets(positive_sample, negtive_sample)#把positive和nagative拼接到一起
        # Generate classifier sample
        print("generating classify sample.....")

        num = len(all_sample)#测量一下获取数据的长度
        start = time.time()#这个是统计计算时间长度写的
        sequence = all_sample["sequence"]#把所有数据的sequence读出来
        peptides = sequence.values.copy().tolist()#把sequience对应的value取出来
        result = all_sample["MIC"]
        type = all_sample["type"]
        output_path = self.generate_example_path + "classify_sample.csv"#输出的文件名称
        cal_pep_des.cal_pep(peptides, sequence, result, type, output_path)#计算结构化数据
        end = time.time()
        print("generate classify feature data cost time:", (end - start) / num)

        print("generating regression sample.....")
        if self.mode == "all":#生成回归用的全部数据
            # generate regression all sample
            self.split_sample(all_sample)
        else:#生成回归用的positive数据
            # generate regression positive sample
            self.split_sample(positive_sample)
        print("regression sample is ok")

    def filter_data_with_str(self, col_name, str, data):#获得列名称中包含"bacterium"（细菌）, "aureus"（金黄色葡萄球菌）的数据

        bool_filter = data[col_name].str.contains(str)
        filter_data = data[bool_filter]
        return filter_data

    def generate_all_peptides(self, data):  # 生成肽类

         data_all = [[], [], []]
         for i in data["sequence"].unique():#分析数据中sequnce列
            data_all[0].append(i)#记录sequence不重复的数值
            log_num = 0
            count = 0
            for i in data[data["sequence"] == i]["value"]:#遍历sequence中每个种类并统计其对应value
                log_num += math.pow(10, i)#计算10**i
                count += 1
            data_all[1].append(float(log_num / count))#log_sum对sequence类内取平均值
            data_all[2].append(1)#似乎是分类标签

         data_all = list(map(list, zip(*data_all)))  # 将列表转置，为后面添加MIC，type的索引准备
         data = pd.DataFrame(data=data_all, columns=["sequence", "MIC", "type"])#在新表里，上述计算为sequence，MIC,type列
         return data

    def data2csv(self, data, file_name):#存储文件

        data.to_csv(file_name, encoding="utf8", index=False)

    def generate_negative_data(self, negative_file_path):#生成negative的数据

        data_negative = pd.read_csv(negative_file_path, encoding="utf8")#读取原始的negative文件
        data_negative = data_negative[~data_negative["Sequence"].str.contains("B|X|Z|O|U")]#只取sequence列不包含"B|X|Z|O|U"的数据
        data_negative.reset_index(drop=True, inplace=True)#重置下标，drop重新设置索引后是否将原索引作为新的一列并入DataFrame,inplace是否在原DataFrame上改动
        data = pd.DataFrame(columns=["sequence", "MIC", "type"])
        for i in range(data_negative.shape[0]):
            data = data.append({"sequence": data_negative["Sequence"][i], "MIC": 8196, "type": 0}, ignore_index=True)

        return data


    def concat_datasets(self, positive_file, negative_file):

        data_concat = pd.concat([positive_file, negative_file], ignore_index=True, axis=0)  # 默认纵向合并0 横向合并1，把sequence,MIC,type分别放在一起，前面带索引
        data_concat = data_concat.sample(frac=1,random_state=None)
        data_concat.reset_index(drop=True, inplace=True)
        return data_concat

    def split_sample(self, sample):
        num = len(sample)
        train_sample = sample[:int(0.8 * num)]
        test_sample = sample[int(0.8 * num):]
        self.data2csv(train_sample, self.generate_example_path + "regression_train_sample.csv")
        self.data2csv(test_sample, self.generate_example_path + "regression_test_sample.csv")
        print('finish')

if __name__ == "__main__":
    grampa_path = "./data/origin_data/grampa.csv"#获取grampa数据集
    negative_file_path = "./data/origin_data/origin_negative.csv"#存储所有认为nagative的样本
    generate_example_path = "./split_files/"#存储划分好的数据
    mode = "all"
    G = GenerateSample(grampa_path, negative_file_path, generate_example_path, mode)
    G()
