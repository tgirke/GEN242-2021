---
title: "Deployment and Maintenance of this Site"
linkTitle: "Installation"
weight: 6
type: docs
---

<br/>
<br/>

> This page provides instructions how to create new deployment instances of this teaching site, and how to configure and customize it. 
> It uses the Docsy theme of the Hugo framework for building content driven websites. 

## Quick start

* Click on the **Use this Template** button.
* Choose a Repository Name
* Click on the **Create repository from template** button.

![](https://raw.githubusercontent.com/dcassol/images/main/usetemplte.gif)

### Usage locally

* Go to your new repository that created from our template `https://github.com/<username>/<repository_name>`
* Click on the **Code** button.
* Copy the URL `https://github.com/<username>/<repository_name>.git`
* Open terminal 

```
git clone --recurse-submodules --depth 1 https://github.com/<username>/<repository_name>.git
cd <repository_name>
```

* Run the website locally

```
hugo server
```

* Run the website locally with `blogdown`

```r
blogdown::serve_site()
```

## Prerequisites and Installation

### Install blogdown and Hugo

#### blogdown

```r
installed.packages("rstudio/blogdown")
# If anything wrong try develop version
remotes::install_github("rstudio/blogdown")
```
#### Hugo

You need a recent extended version (we recommend version 0.79.0 or later) of Hugo 
to do local builds and previews of sites that use Docsy.

It is recommended to install `Hugo` from R for working with {blogdown}

```r
blogdown::install_hugo(extended = TRUE)
```
or from commandline

```bash
wget https://github.com/gohugoio/hugo/releases/download/v0.79.0/hugo_extended_0.79.0_Linux-64bit.deb
sudo dpkg -i  hugo_extended_0.79.0_Linux-64bit.deb 
hugo version
```

For `Windows` and `macOS` please see instructions [here](https://www.docsy.dev/docs/getting-started/). 

#### Install PostCSS


To build or update your site's CSS resources, you also need `PostCSS` to create the final assets. If you need to install it, you must have a recent version of `NodeJS` installed on your machine so you can use `npm`, the Node package manager. By default `npm` installs tools under the directory where you run `npm` install:

```bash
sudo npm install -D autoprefixer
sudo npm install -D postcss-cli
# Starting in version 8 of postcss-cli, you must also separately install postcss:

sudo npm install -D postcss

# go to your website directory
cd <repository_name>
npm audit fix
```

### Run the website locally 

#### with {blogdown}

* Open R in console or Rstudio

This repo contains a `.Rprofile` file that will automatically serve the site
for you R starting directory is this newly cloned repository. Otherwise, 
change working directory to the repository and run:

```r
blogdown::serve_site()
```

You should see a website is opened in your local browser or Rstudio viewer.

#### Run the website locally on the terminal

```bash
cd YOUR_NEW_REPO_PATH
hugo server
```
