from openquake.baselib.general import git_suffix

__version__ = '1.3'
__version__ += git_suffix(__file__)

__import__('pkg_resources').declare_namespace(__name__)
