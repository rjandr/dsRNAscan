import os
import sys
import platform
from setuptools import setup, find_packages

# Version information
__version__ = '0.5.0'

def get_platform_binary():
    """Determine which pre-compiled einverted binary to use"""
    system = platform.system().lower()
    machine = platform.machine().lower()
    
    if system == 'linux':
        if machine in ['x86_64', 'amd64']:
            return 'einverted_linux_x86_64'
        elif machine in ['aarch64', 'arm64']:
            return 'einverted_linux_arm64'
        elif machine in ['i386', 'i686']:
            return 'einverted_linux_i386'
    elif system == 'darwin':  # macOS
        if machine in ['arm64', 'aarch64']:
            return 'einverted_macos_arm64'
        elif machine in ['x86_64', 'amd64']:
            return 'einverted_macos_x86_64'
    elif system == 'windows':
        return 'einverted_windows_x86_64.exe'
    
    # Fallback to generic Linux
    return 'einverted_linux_x86_64'

def setup_einverted():
    """Set up the appropriate einverted binary"""
    setup_dir = os.path.dirname(os.path.abspath(__file__))
    tools_dir = os.path.join(setup_dir, 'dsrnascan', 'tools')
    platform_binaries_dir = os.path.join(tools_dir, 'platform_binaries')
    
    # Determine which binary to use
    binary_name = get_platform_binary()
    source_binary = os.path.join(platform_binaries_dir, binary_name)
    
    # Target location (without platform suffix)
    if binary_name.endswith('.exe'):
        target_binary = os.path.join(tools_dir, 'einverted.exe')
    else:
        target_binary = os.path.join(tools_dir, 'einverted')
    
    # Check if source binary exists
    if os.path.exists(source_binary):
        print(f"Using pre-compiled einverted: {binary_name}")
        # Copy to target location if needed
        if not os.path.exists(target_binary) or os.path.getmtime(source_binary) > os.path.getmtime(target_binary):
            import shutil
            shutil.copy2(source_binary, target_binary)
            # Make executable on Unix-like systems
            if not binary_name.endswith('.exe'):
                os.chmod(target_binary, 0o755)
            print(f"✓ Installed einverted binary for {platform.system()} {platform.machine()}")
    else:
        print(f"WARNING: No pre-compiled binary found for {platform.system()} {platform.machine()}")
        print(f"  Looking for: {source_binary}")
        print("  einverted functionality may not be available")
        
        # List available binaries
        if os.path.exists(platform_binaries_dir):
            available = os.listdir(platform_binaries_dir)
            if available:
                print(f"  Available binaries: {', '.join(available)}")

# Set up einverted when installing
if any(cmd in sys.argv for cmd in ['install', 'build', 'bdist_wheel', 'develop']):
    setup_einverted()

# Read long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='dsrnascan',
    version=__version__,
    author='Bass Lab',
    author_email='',
    description='A tool for genome-wide prediction of double-stranded RNA structures',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/Bass-Lab/dsRNAscan',
    project_urls={
        "Bug Tracker": "https://github.com/Bass-Lab/dsRNAscan/issues",
        "Documentation": "https://github.com/Bass-Lab/dsRNAscan/blob/main/README.md",
        "Source Code": "https://github.com/Bass-Lab/dsRNAscan",
    },
    packages=['dsrnascan'],
    package_data={
        'dsrnascan': [
            'tools/*',
            'tools/platform_binaries/*',
            'models/stability_rf_depth5/*',
            'models/probing_rf/*',
        ],
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.8',
    install_requires=[
        'biopython>=1.78',
        'numpy>=1.19',
        'pandas>=1.1',
        'ViennaRNA>=2.4',
        'psutil>=5.8',
        'tqdm>=4.0',
        'ydf>=0.9.0',
        'scikit-learn>=1.0',
    ],
    extras_require={
        'mpi': ['mpi4py>=3.0', 'parasail>=1.2'],
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=22.0',
            'flake8>=4.0',
            'mypy>=0.900',
        ],
    },
    entry_points={
        'console_scripts': [
            'dsrnascan=dsrnascan:main',
            'dsrna-browse=dsrnascan.browse_results:main',
        ],
    },
)