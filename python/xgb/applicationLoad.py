import sys
import pandas
import cPickle as pickle
import itertools
import collections
import traceback
import xgboost


def load(pklfile):
  f = None
  pkldata = None
  try:
    f = open(pklfile,'rb')
  except IOError as e:
    print('Couldnt open or read from file (%s).' % e)
  else:
      try:
        pkldata = pickle.load(f)
      except pickle.UnpicklingError as e: # normal, somewhat expected
        try:
          model = pkldata.booster().get_dump() # this only tests load was ok
        except (AttributeError,  EOFError, ImportError, IndexError) as e:
          print(traceback.format_exc(e))
        except Exception as e:
          print(traceback.format_exc(e))
      f.close()
  return pkldata


def evaluate(vec_values, vec_names, pkldata):
  new_dict = collections.OrderedDict(itertools.izip(vec_names, vec_values))
  data = pandas.DataFrame(columns = list(new_dict.keys()))
  data = data.append(new_dict, ignore_index = True)
  if 'f0' in pkldata.get_booster().feature_names: # For nameless features
    data_to_use = data[data.columns.values.tolist()].values
  else: # For named features
    data_to_use = data
  result = -20
  if 'XGBClassifier' in str(type(pkldata)):
    try:
      proba = pkldata.predict_proba(data_to_use)
    except:
      print('Caught error:', sys.exc_info())
    else:
      result = proba[:,1][0]
  else:
    try:
      matrix = xgboost.DMatrix(data,feature_names=new_dict.keys())
      proba = pkldata.predict(matrix)
    except:
      print('Caught error:', sys.exc_info())
    else:
      result = proba[:,1][0]
  return result
