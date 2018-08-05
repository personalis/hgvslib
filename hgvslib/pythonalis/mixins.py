import os
import sys

import inspect

# TODO standardize on "Mixin" suffix or no suffix

class Registry(object):
   __registries = {}

   @classmethod
   def _add_registry(cls, registry_name, registry_obj):
      if registry_name in cls.__registries:
         sys.stderr.write('Overwriting registry "{0}"\n'.format(registry_name))
      cls.__registries[registry_name] = registry_obj

   @classmethod
   def _has_registry(cls, registry_name):
      return registry_name in cls.__registries

   @classmethod
   def registry_lookup(cls, registry_name, key):
      if not cls._has_registry(registry_name):
         sys.stderr.write('Registry "{0}" not found.\n'.format(registry_name))

      elif not cls.registry_exists(registry_name, key):
         sys.stderr.write('Registry "{0}" does not have key "{1}"\n'.format(registry_name, key))

      return cls.__registries.get(registry_name, {}).get(key, None)

   @classmethod
   def registry_exists(cls, registry_name, key):
      if not cls._has_registry(registry_name):
         sys.stderr.write('Registry "{0}" not found.\n'.format(registry_name))

      return key in cls.__registries.get(registry_name, {})

   @classmethod
   def registry_insert(cls, registry_name, key, value):
      cls.__registries[registry_name][key] = value

   @classmethod
   def registry_update(cls, registry_name, val_dict):
      if not cls._has_registry(registry_name):
         sys.stderr.write('Registry "{0}" not found.\n'.format(registry_name))
         return False

      cls.__registries[registry_name].update(val_dict)
      return True