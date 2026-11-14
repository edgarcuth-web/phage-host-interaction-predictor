"""Independently check saved Kaneko features, scores, metrics and intervals.

No fitting functions or project scoring/metric helpers are imported.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
from itertools import product
import json
from math import comb, exp, isfinite
from pathlib import Path
import re

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = "kaneko-2025:complete-spot-clearing-grade3:v1"
SECONDARY = "kaneko-2025:visible-spot-clearing-grade-ge2:v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def close(left, right, name):
    if left is None or right is None:
        require(left is None and right is None, "Undefined-value mismatch: "+name)
    else:
        require(isfinite(float(left)) and isfinite(float(right)) and abs(float(left)-float(right)) <= 1e-12,
                "Numerical mismatch: "+name)


def verify_binding(path, expected):
    data=Path(path).read_bytes()
    require(len(data)==expected["size_bytes"] and sha(data)==expected["sha256"], "Changed bytes: "+str(path))


def source_fourmers(payload, words):
    """Use direct forward/reverse 2-bit indices rather than the project bin map."""
    contigs=[re.sub(r"\s+", "", record.split("\n",1)[1]).upper()
             for record in payload.decode("ascii").split(">")[1:]]
    require(contigs and all(contigs), "Empty FASTA source")
    lookup=np.full(256,4,dtype=np.uint16)
    for value,letter in enumerate(b"ACGT"):
        lookup[letter]=value
    counts=np.zeros(256,dtype=np.int64)
    for sequence in contigs:
        encoded=lookup[np.frombuffer(sequence.encode(),dtype=np.uint8)]
        if len(encoded)<4:
            continue
        a,b,c,d=encoded[:-3],encoded[1:-2],encoded[2:-1],encoded[3:]
        valid=(a<4)&(b<4)&(c<4)&(d<4)
        forward=64*a+16*b+4*c+d
        reverse=64*(3-d)+16*(3-c)+4*(3-b)+(3-a)
        counts+=np.bincount(np.minimum(forward,reverse)[valid],minlength=256)
    word_indices=[sum("ACGT".index(base)*4**(3-i) for i,base in enumerate(word)) for word in words]
    require(counts.sum()>0 and counts[word_indices].sum()==counts.sum(), "Invalid independent k-mer counts")
    return counts[word_indices]/counts.sum()


def independent_metrics(labels, scores):
    blocks=[]
    for value in sorted(set(scores),reverse=True):
        selected=[int(y) for y,s in zip(labels,scores) if s==value]
        blocks.append((len(selected),sum(selected)))
    positives=sum(labels)
    hits=[0.0]
    recalls={}
    for k in range(1,len(scores)+1):
        remaining=k
        no_hit=1.0
        recovered=0.0
        for n,p in blocks:
            take=min(n,remaining)
            remaining-=take
            no_hit*=comb(n-p,take)/comb(n,take) if take<=n-p else 0.0
            recovered+=take*p/n
        hits.append(1-no_hit)
        if k<=5:
            recalls[k]=recovered/positives if positives else None
    result={**{f"hit@{k}":hits[k] for k in range(1,6)},
            **{f"recall@{k}":recalls[k] for k in range(1,6)},
            "mrr":sum((hits[k]-hits[k-1])/k for k in range(1,len(scores)+1)),
            "capped_tests@2":1+sum(1-hits[k] for k in range(1,2)),
            "capped_tests@5":1+sum(1-hits[k] for k in range(1,5))}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit-dir",type=Path,required=True)
    parser.add_argument("--evaluation-dir",type=Path,required=True)
    parser.add_argument("--sequence-root",type=Path,required=True)
    parser.add_argument("--grades",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    require(not args.out.exists(),"Verification output must be new")
    protocol=json.loads((ROOT/"data/studies/results/kaneko_baseline_freeze.v1/protocol.json").read_text())
    design=json.loads((ROOT/protocol["design_path"]).read_text())
    library=protocol["library"]
    host_rows=rows(ROOT/"data/studies/results/kaneko_baseline_freeze.v1/hosts.csv")
    host_info={r["host_id"]:r for r in host_rows}
    fit=json.loads((args.fit_dir/"FIT.json").read_text())
    report=json.loads((args.evaluation_dir/"report.json").read_text())
    for name,item in fit["files"].items():
        verify_binding(args.fit_dir/name,item)
    for name,item in json.loads((args.evaluation_dir/"OUTPUTS.json").read_text()).items():
        verify_binding(args.evaluation_dir/name,item)
    require(sha(args.grades.read_bytes())==design["source_grades"]["sha256"],"Source grades changed")
    grade_map={(r["host_label"],r["phage_label"]):int(r["source_grade"]) for r in rows(args.grades)}
    words=sorted({min(w,w.translate(str.maketrans("ACGT","TGCA"))[::-1])
                  for w in map("".join,product("ACGT",repeat=4))})
    feature_rows=rows(args.fit_dir/"host_features.csv")
    features={r["host_id"]:np.array([float(r[w]) for w in words]) for r in feature_rows}
    require(set(features)==set(host_info) and len(features)==143,"Feature host identities differ")
    source=json.loads((ROOT/design["inputs"]["sequence_intake"]["path"]).read_text())
    for item in source["host_assemblies"]:
        path=args.sequence_root/item["fasta"]["filename"]
        verify_binding(path,item["fasta"])
        vector=source_fourmers(gzip.decompress(path.read_bytes()),words)
        require(np.array_equal(vector,features[item["source_label"]]),"Source feature vector differs")
    model=json.loads((args.fit_dir/"host_4mer_model.json").read_text())
    prevalence=json.loads((args.fit_dir/"prevalence_model.json").read_text())
    train=sorted(h for h,r in host_info.items() if r["partition"]=="train")
    require(model["training_hosts"]==prevalence["training_hosts"]==train and len(train)==123,"Training hosts differ")
    require(set(model["phages"])==set(prevalence["phages"])==set(library),"Model library differs")
    require(model["endpoint_id"]==prevalence["endpoint_id"]==PRIMARY,"Model endpoint differs")
    require(model["score_semantics"]==prevalence["score_semantics"]=="uncalibrated_source_grade3_complete_spot_clearing",
            "Historical score descriptor propagated")
    matrix=np.array([features[h] for h in train])
    means,scales=matrix.mean(axis=0),matrix.std(axis=0)
    scales[scales==0]=1
    np.testing.assert_allclose(means,model["preprocessing"]["means"],rtol=0,atol=1e-14)
    np.testing.assert_allclose(scales,model["preprocessing"]["scales"],rtol=0,atol=1e-14)
    means=np.array(model["preprocessing"]["means"])
    scales=np.array(model["preprocessing"]["scales"])
    for phage in library:
        positives=sum(grade_map[h,phage]==3 for h in train)
        expected=(positives+1)/(len(train)+2)
        close(prevalence["phages"][phage]["score"],expected,"training prevalence")
        if model["phages"][phage]["kind"]=="single_class_prevalence":
            close(model["phages"][phage]["score"],expected,"single-class fallback")
        require(model["phages"][phage]["training_rows"]==123,"Per-phage fit scope differs")
    score_rows=rows(args.evaluation_dir/"scores.csv")
    require(len(score_rows)==1680,"Score row count differs")
    score_map={(r["partition"],r["method"],r["host_id"],r["phage_id"]):float(r["score"]) for r in score_rows}
    require(len(score_map)==1680 and all(r["score_endpoint_id"]==PRIMARY for r in score_rows),"Duplicate/wrong-endpoint score")
    maximum_score_error=0.0
    with threadpool_limits(limits=1):
        for key,value in score_map.items():
            part,method,host,phage=key
            require(host_info[host]["partition"]==part and part in {"validation","test"},"Score host leakage")
            if method=="prevalence":
                expected=prevalence["phages"][phage]["score"]
            else:
                require(method=="host_4mer","Unexpected method")
                local=model["phages"][phage]
                if local["kind"]=="single_class_prevalence":
                    expected=local["score"]
                else:
                    z=float(((features[host]-means)/scales)@np.array(local["coefficients"])+local["intercept"])
                    expected=1/(1+exp(-z)) if z>=0 else exp(z)/(1+exp(z))
            maximum_score_error=max(maximum_score_error,abs(value-expected))
            close(value,expected,"saved candidate score")
    independent={}
    per_host=rows(args.evaluation_dir/"per_host_metrics.csv")
    metric_checks=0
    for row in per_host:
        endpoint,part,method,host=(row[k] for k in ("evaluation_endpoint_id","partition","method","bacterium"))
        minimum=3 if endpoint==PRIMARY else 2
        require(endpoint in {PRIMARY,SECONDARY} and row["score_endpoint_id"]==PRIMARY,"Unexpected metric endpoint")
        target=[int(grade_map[host,p]>=minimum) for p in library]
        scores=[score_map[part,method,host,p] for p in library]
        result=independent_metrics(target,scores)
        result["eligible"]=bool(sum(target))
        require(int(row["positive_candidates"])==sum(target),"Positive support differs")
        require(row["host_group"]==host_info[host]["host_group"],"Metric host group differs")
        for metric,value in result.items():
            if metric=="eligible":
                require(row["eligible"]==str(value),"Eligibility differs")
            else:
                saved=None if row[metric]=="" else float(row[metric])
                close(saved,value,"per-host "+metric)
                metric_checks+=1
        independent[endpoint,part,method,host]=result
    require(len(per_host)==len(independent)==80,"Per-host result scope differs")
    contrast_rows=json.loads((args.evaluation_dir/"contrasts.json").read_text())
    for contrast in contrast_rows:
        endpoint,part,metric=(contrast[k] for k in ("evaluation_endpoint_id","partition","metric"))
        hosts=sorted(h for h,r in host_info.items() if r["partition"]==part)
        names=sorted({host_info[h]["host_group"] for h in hosts})
        group_sums,group_counts=[],[]
        for group in names:
            members=[h for h in hosts if host_info[h]["host_group"]==group]
            values=[]
            for h in members:
                a,b=independent[endpoint,part,"host_4mer",h],independent[endpoint,part,"prevalence",h]
                if not metric.startswith("recall@") or a["eligible"]:
                    values.append(a[metric]-b[metric])
            group_sums.append(sum(values))
            group_counts.append(len(values))
        sums,counts=np.array(group_sums),np.array(group_counts)
        draws=np.random.default_rng(20261114).integers(len(names),size=(2000,len(names)))
        denominators=counts[draws].sum(axis=1)
        valid=denominators>0
        values=sums[draws].sum(axis=1)[valid]/denominators[valid]
        interval=np.quantile(values,[0.025,0.975]).tolist() if len(values) else None
        close(contrast["delta"],float(sums.sum()/counts.sum()) if counts.sum() else None,"contrast point")
        require(contrast["valid_resamples"]==int(valid.sum()) and contrast["undefined_resamples"]==int((~valid).sum()),
                "Bootstrap support differs")
        if interval is None:
            require(contrast["interval"] is None,"Undefined interval differs")
        else:
            for a,b in zip(contrast["interval"],interval):
                close(a,b,"group bootstrap interval")
    require(len(contrast_rows)==44,"Contrast scope differs")
    phage_rows=rows(args.evaluation_dir/"per_phage_diagnostics.csv")
    expected_selection=defaultdict(lambda: [0.0,0.0])
    for endpoint,minimum in [(PRIMARY,3),(SECONDARY,2)]:
        for part in ["validation","test"]:
            for method in ["prevalence","host_4mer"]:
                for host in sorted(h for h,r in host_info.items() if r["partition"]==part):
                    values=[score_map[part,method,host,p] for p in library]
                    remaining=2
                    for value in sorted(set(values),reverse=True):
                        positions=[i for i,x in enumerate(values) if x==value]
                        mass=min(remaining,len(positions))/len(positions)
                        remaining-=min(remaining,len(positions))
                        for i in positions:
                            item=expected_selection[endpoint,part,method,library[i]]
                            item[0]+=mass
                            item[1]+=mass*int(grade_map[host,library[i]]>=minimum)
    for row in phage_rows:
        key=tuple(row[k] for k in ("evaluation_endpoint_id","partition","method","phage_id"))
        for saved,expected in zip([row["expected_top2_selections"],row["expected_positive_top2_selections"]],expected_selection[key]):
            close(float(saved),expected,"per-phage tie selection")
    require(len(phage_rows)==336,"Per-phage diagnostic scope differs")
    hit=report["primary_hit_contrast"]["interval"]
    cost=report["primary_cost_contrast"]["interval"]
    require(report["decisions"]=={"hit_gain_supported":hit[0]>0,"cost_reduction_supported":cost[1]<0,
             "joint_benefit_supported":hit[0]>0 and cost[1]<0},"Acceptance decision differs")
    require(fit["held_out_score_runs"]==0 and report["fit_calls_during_evaluation"]==0,"Phase boundary differs")
    result={"schema_version":1,"study_id":design["study_id"],"status":"independent_numeric_checks_passed",
       "fit_commit":report["frozen_fit_commit"],"implementation_commit":fit["implementation_commit"],
       "source_host_feature_vectors_verified":143,"feature_values_verified":143*136,
       "source_models_verified":42,"candidate_score_values_verified":len(score_rows),
       "maximum_score_absolute_error":maximum_score_error,"per_host_metric_values_verified":metric_checks,
       "bootstrap_intervals_verified":len(contrast_rows),"per_phage_rows_verified":len(phage_rows),
       "refits_performed_by_verifier":0,"decisions":report["decisions"],
       "verifier_sha256":sha(Path(__file__).read_bytes()),"independent_release_validation_established":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open("x") as stream:
        stream.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result),flush=True)


if __name__=="__main__":
    main()
