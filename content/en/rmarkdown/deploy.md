---
title: "Deploy"
linkTitle: "Deploy"
weight: 4
type: docs
---

### Deploy to GitHub Pages

Key feature: Fully automatic, **no need to build website locally**. Yes, you can skip `blogdown::build_site()` step locally.

#### 1. Add all R and system packages to `/deps.yaml` file that are required to render your `Rmarkdown` files

* Specify any CRAN, Bioconductor or Github packages you have used, in `yaml` array format. Github packages need to be `user_name/repo_name`.
* The site building is happening on Github by an Ubuntu system, so specify any system packages to install that are required by your R packages. Only a single line, use space to separate each system dependencies.

#### 2. Change base URL in config file

Change your `baseURL` in the `/config.toml`. 

* If this is your first github website, change the baseURL = "/" to your github base `url`, like "`https://USERNAME.github.io/`".
* If not, baseURL = "/" to your repo's url, like "`https://USERNAME.github.io/REPO_NAME/`".
* If you use custom domain, baseURL = "/" to your custom url, like "`https://MY_DOMAIN_NAME.com/`".

#### 3. Create a new branch named `gh-pages`

```bash
git checkout -b gh-pages
rm -rv !(.git|README.md)
git add . && git commit -m "gh-pages" && git push --set-upstream origin gh-pages
git checkout main
git add . && git commit -m "pg_build" && git push
git submodule update --init --recursive
```
Or go to Github website create a branch by [clicking](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-and-deleting-branches-within-your-repository).

You do *not* need to push any content to this `gh-pages` branch. This branch is automatically managed by Github Actions. 

#### 4. Github Setting 

* Click on **Setting**
* click on **GitHub Pages** 
* Select Branch *gh-pages* and folder */(root)* and click on **Save**
* Click on **Enforce HTTPS** 

#### 5. Trigger page build

* To trigger a build and deployment event, in the `git commit` message, include the word "**pg_build**", 
e.g. `git commit -m "blabla pg_build"`. This will run the full rendering process (rendering and deploy the website may still take a few minutes).
* If you have build the site locally and updated the `/public` locally, in your commit message, include the word "**no_render**" will skip the site rendering + build process. Github Actions will directly use `/public` to host the website.

#### Details

* The first deployment can take some time because Github Actions need to install these packages. 
After first successful deployment with Github Actions, all packages will be cached and 
no need to install again later. So later deployment will be much faster.

 
## Deploy to Netlify

[![Netlify Status](https://api.netlify.com/api/v1/badges/2de20eae-a002-40ef-8c96-1fe54f9528d4/deploy-status)](https://app.netlify.com/sites/docsy-rmarkdown/deploys)

**_Important:_** :

1. Before running `blogdown::build_site()`, in `/config.toml` change `baseURL` to "/" `baseURL="/"`
1. Before deployment, in `/netlify.toml` change `publish` to "public" `publish = "public"`

### Deploy to Netlify with GitHub

* Go online to [Netlify.com](https://www.netlify.com/).
* Click on the **Sign Up** button.
* It is recommended to sign up using your existing GitHub account, so select **GitHub**, and click to **Authorize Netlify**.
* Click on the **New site from Git** button.
* Click on the **GitHub** button.
* Select the repository you just created it.
* Click on **Show advanced** and then **New Variable**:
  * Specify `HUGO_VERSION` as the **Key** for the new variable, and `0.79.0` for the **Value**.
* Click on **Deploy site**.

### Deploy to Netlify button

[![Deploy to Netlify](https://www.netlify.com/img/deploy/button.svg)](https://app.netlify.com/start/deploy?repository=https://github.com/dcassol/docsy_Rmarkdown)

> **_NOTE:_** Gitib does not allow third party tools to create/modify Github actions workflows and this template contains workflows, so 
clicking the button above **will not copy code** to your new repository, but it will publish website to netlify. Follow the instructions above to 
create a repository with code. 

Click on the button above and follow these instructions:

* Click **Connect to GitHub**
* Enter **Login** and **Password**
* Choose Repository Name and **Save & Deploy**

<details>
<summary><b>
Add your own status badge
</b></summary

You can get a status badge for you website by going to Netlify website, login, go to `Site settings` > `Status badges`, copy the status badge. 

</details>
