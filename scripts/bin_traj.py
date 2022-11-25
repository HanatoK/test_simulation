#!/usr/bin/env python3
import pandas as pd
from histogram.build_histogram_from_traj import BuildHistogramFromTraj
from histogram.reweight import convert_probability_to_pmf


def main():
    data = pd.read_csv('../build/XYZ.traj', delimiter=r'\s+', comment='#', header=None)
    build_histogram = BuildHistogramFromTraj(json_file='axis.json', position_cols=[1, 2])
    build_histogram.read_pandas(data)
    with open('hist.dat', 'w') as f_output:
        convert_probability_to_pmf(build_histogram.get_histogram(), 300.0*0.0019872041).write_to_stream(stream=f_output)


if __name__ == '__main__':
    main()
