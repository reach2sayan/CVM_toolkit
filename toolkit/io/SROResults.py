from __future__ import annotations
from dataclasses import dataclass
from typing import Union
import numpy as np
import pandas as pd

@dataclass(kw_only=True, eq=False, order=False)
class SROResults:

    phase: str
    structure: str
    norm_constrained: bool = False
    _result: pd.DataFrame = pd.DataFrame(columns=['phase','structure','temperature',
                                                  'F_sqs','F_opt','F_ord','F_rnd',
                                                  'constr_tol','correlations'
                                                 ]
                                        )

    @property
    def result(self: SROResults) -> pd.DataFrame:
        return self._result

    @result.setter
    def result(self: SROResults,
               df: Union[str,pd.DataFrame,list[dict]],
              ) -> pd.DataFrame:
        if isinstance(df, str):
            self._result = pd.read_csv(df)
        elif isinstance(df, pd.DataFrame):
            self._result = df.copy()
        elif isinstance(df, list):
            if all(isinstance(x, dict) for x in df):
                self._result = pd.DataFrame.from_dict(df,orient='columns')
        else:
            raise ValueError('result should either be a string reference to a csv file or a pandas dataframe or a list of dictionaries.\
                             The columns of the csv file or the keys of the dicionaries, should be\n\
                             "phase","structure","temperature","F_sqs",F_opt","F_ord","F_rnd","constr_tol","correlations"'
                            )

    def save_to_file(self: SROResults, outfile) -> None:
        self._result.to_csv(outfile,index=False)
