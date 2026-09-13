from pathlib import Path
import json
import tempfile
import numpy as np
from solve import integrate, read_inputs
from export_result2_full import export_run, audit_book


def main():
    arrays,_=read_inputs()
    rows=[]
    for q in (1,3,4):
        run=integrate(q,80,arrays,horizon=600,stop_at_threshold=False)
        sparse=run.summary()
        dense=run.summary(sample_step_s=1,chunk_size=37)
        assert sparse['diagnostic_sampling']['time_count']==250
        assert dense['diagnostic_sampling']['time_count']==601
        assert not dense['diagnostic_sampling']['continuous_extrema_claimed']
        assert not {'T_max','T_min','C_max','C_min','moisture_balance_max','max_above_center'} & dense.keys()
        y=run.values(np.arange(601));m=81
        expected={'T_max_sampled':y[:m].max(),'T_min_sampled':y[:m].min(),
                  'C_max_sampled':y[m:2*m].max(),'C_min_sampled':y[m:2*m].min(),
                  'moisture_balance_max_sampled':np.max(np.abs(2*run.model.grid.w@y[m:2*m]+y[-1]-2.55))}
        for key,value in expected.items():
            assert abs(dense[key]-value)<1e-12,(q,key,dense[key],value)
        assert dense['moisture_balance_max_sampled']<1e-7
        for bad in (0,-1,float('nan'),float('inf')):
            try:run.summary(sample_step_s=bad)
            except ValueError:pass
            else:raise AssertionError('Invalid step accepted')
        rows.append(dict(question=q,diagnostics=dense))
        if q==3:
            out=Path(tempfile.mkdtemp(prefix='smoke_output_',dir=Path(__file__).resolve().parent))
            export_run(run,out,chunk=137)
            audit=audit_book(out/'result2.xlsx',600)
            assert len(audit)==2 and all(a['data_rows']==600 for a in audit)
    print(json.dumps(dict(status='PASS',runs=rows),ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
