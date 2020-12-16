conda activate bio

cd ~/python_work/nCov/bin

python3 Iteration.py -q ../data/PangolinGD_ORF1ab.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.ORF1ab.result

python3 Iteration.py -q ../data/PangolinGD_S.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.S.result

python3 Iteration.py -q ../data/PangolinGD_E.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.E.result

python3 Iteration.py -q ../data/PangolinGD_M.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.M.result

python3 Iteration.py -q ../data/PangolinGD_N.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.N.result

python3 Iteration.py -q ../data/PangolinGX_ORF1ab.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.ORF1ab.result

python3 Iteration.py -q ../data/PangolinGX_S.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.S.result

python3 Iteration.py -q ../data/PangolinGX_E.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.E.result

python3 Iteration.py -q ../data/PangolinGX_M.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.M.result

python3 Iteration.py -q ../data/PangolinGX_N.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.N.result

