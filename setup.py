from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='crysfipy',
      version='0.5.2',
      description='Crystal field suite for python.',
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
      ],
      keywords=['physics', 'crystal field', 'magnetism'],
      url='https://bitbucket.org/cermak/crysfipy',
      download_url='https://bitbucket.org/cermak/crysfipy/get/0.5.tar.gz',
      author='Petr Čermák',
      author_email='pcermak@live.com',
      license='GPLv3',
      packages=['crysfipy'],
      install_requires=[
          'numpy',
      ],
      include_package_data=True,
      zip_safe=False)