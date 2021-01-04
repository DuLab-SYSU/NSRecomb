conda activate bio

cd ~/python_work/nCov/bin

python3 Iteration.py -q ../data/PangolinGD.ORF1ab.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.ORF1ab.result

python3 Iteration.py -q ../data/PangolinGD.S.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.S.result

python3 Iteration.py -q ../data/PangolinGD.E.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.E.result

python3 Iteration.py -q ../data/PangolinGD.M.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.M.result

python3 Iteration.py -q ../data/PangolinGD.N.fasta -b ../data/PangolinGD.backbone  -n ../data/PangolinGD.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGD.N.result

python3 Iteration.py -q ../data/PangolinGX.ORF1ab.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.ORF1ab.result

python3 Iteration.py -q ../data/PangolinGX.S.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.S.result

python3 Iteration.py -q ../data/PangolinGX.E.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.E.result

python3 Iteration.py -q ../data/PangolinGX.M.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.M.result

python3 Iteration.py -q ../data/PangolinGX.N.fasta -b ../data/PangolinGX.backbone  -n ../data/PangolinGX.acc -db ../db/all-cov -ws 501 -nc 10 -o ../result/PangolinGX.N.result

