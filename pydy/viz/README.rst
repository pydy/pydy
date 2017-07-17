This is how you setup a development install of an nbextension for a conda env::

   $ jupyter nbextension install <directory-with-nbextension-js-files> --sys-prefix --symlink --destination=<new-directory-name>
   $ jupyter nbextension enable <new-directory-name>/<name-of-js-file-without-extension> --sys-prefix

For example::

   $ ls pydy/viz/nbextension/
   pydyviz.js
   $ jupyter nbextension install pydy/viz/nbextension/ --sys-prefix --symlink --destination="pydyviz"
   $ jupyter nbextension enable pydyviz/pydyviz --sys-prefix
