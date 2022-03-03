import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
        name='CAPRToolbox',
        version = '0.0.1',
        author = 'Yang Wang',
        author_email = 'yplsw90@outlook.com',
        description = 'Testing installation of Package',
        long_description = long_description,
        long_description_content_type = "text/markdown",
        url = 'http://github.com/YangChemE/CAPRToolBox',
        project_urls = {
            "Bug Tracker": "https://github.com/YangChemE/CAPRToolbox/issues"
        },
        license='MIT',
        packages=['CAPRToolbox'],
        install_requires=['request'],
)
