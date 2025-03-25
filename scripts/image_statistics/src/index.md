---
title: Image Statistics
toc: false
---

<!-- markdownlint-disable -->


# Image Statistics

<!-- Set up useful functions -->

```js
import {DateTime} from "npm:luxon";
```

```js
function formatTimestamp(datetime) {
  if(!datetime) return null;
  return DateTime.fromMillis(datetime).toLocaleString(DateTime.DATETIME_SHORT)
}
```

<!-- Load the data -->

```js
const db = DuckDBClient.of({
  logs: FileAttachment("./data/logs.parquet"),
  images: FileAttachment("./data/images.parquet")
});
```



<!-- First table view, just a list of repositories -->

```js
const repositoryList = await db.sql`SELECT distinct repository from images`;
const repositoryTable = Inputs.table(repositoryList, {
  format: {
    repository: (repo) => html`<a href=https://console.cloud.google.com/artifacts/docker/cpg-common/australia-southeast1/${repo}?project=cpg-common target=_blank>${repo}</a>`
  }
});
const repositoryTableView = view(repositoryTable);
```

<!-- Image search view, queries a list of images and some associated stats and populates the search box -->

```js
const selected = repositoryTableView.map(rr => `'${rr.repository}'`).join(', ');

const imagesInSelectedRegistries = await db.query(`
  select
    i.short_path as image,
    count(distinct i.digest) filter (where i.status = 'active') as versions,
    count(distinct i.digest) filter (where i.status = 'archived') as archived_versions,
    count(l.timestamp) as total_pulls,
    count(l.timestamp) filter (where i.status = 'active' and l.timestamp > (current_date - interval 1 month)) as pulls_in_last_month,
    max(l.timestamp) as most_recent_pull,
    min(l.timestamp) as first_pull

  from images i
  left join logs l
  on l.full_path = i.full_path
  where true
  ${selected && `and i.repository IN (${selected})`}

  group by 1
  order by 1
`);

const imageSearch = Inputs.search(imagesInSelectedRegistries);
const imageSearchView = view(imageSearch);
```


<!--
  Image list table
  This needs to be a separate block to the search view due to how reactivity works in observable
-->

```js
const imageList = Inputs.table(
  imageSearchView,
  {
    layout: 'auto',
    format: {
      most_recent_pull: (dd) => formatTimestamp(dd),
      first_pull: (dd) => formatTimestamp(dd),
      image: (image) => html`<a href=https://console.cloud.google.com/artifacts/docker/cpg-common/australia-southeast1/${image}?project=cpg-common target=_blank>${image}</a>`

    }
  }
);
const imageListView = view(imageList);

```


<!-- Allow selecting the status of images -->
```js
const statusSelect = Inputs.radio(["all", "active", "archived"], {label: "status", value: "all"});
const statusSelectView = view(statusSelect);
```

<!-- Filter images by status and populate image version search -->
```js
const selectedImages = imageListView.map(rr => `'${rr.image}'`).join(', ');
const statusSelectQuery = statusSelectView === 'all' ? '' : `and i.status = '${statusSelectView}'`;


const selectedImageVersions = await db.query(`
  select
    i.short_path as image,
    i.tags,
    i.digest,
    i.status,
    i.size_bytes / 1024 / 1024 / 1024 as size_gb,
    i.upload_time,
    i.build_time,
    i.update_time,
    count(l.timestamp) as total_pulls,
    count(l.timestamp) filter (where i.status = 'active' and l.timestamp > (current_date - interval 1 month)) as pulls_in_last_month,
    max(l.timestamp) as most_recent_pull
  from images i
  left join logs l
  on l.full_path = i.full_path
  where TRUE
  ${selectedImages && `AND repository || '/' || name IN (${selectedImages})`}
  ${selected && `AND repository IN (${selected})`}
  ${statusSelectQuery}
  group by all
  order by i.short_path, upload_time desc
`);

const imageVersionSearch = Inputs.search(selectedImageVersions);
const imageVersionSearchView = view(imageVersionSearch);

```

<!-- Image version table -->
```js


const imageVersions = Inputs.table(imageVersionSearchView, {
  layout: 'auto',
  format: {
    most_recent_pull: (dd) => formatTimestamp(dd),
    build_time: (dd) => formatTimestamp(dd),
    upload_time: (dd) => formatTimestamp(dd),
    update_time: (dd) => formatTimestamp(dd),
    image: (image) => html`<a href=https://console.cloud.google.com/artifacts/docker/cpg-common/australia-southeast1/${image}?project=cpg-common target=_blank>${image}</a>`,
    digest: (digest, index, items) => {
      // Transform the repo name to add back in the `-archive` suffix so that links work as expected.
      const row = items[index];
      const [repo, ...imageNameList] = row.image.split('/');
      const image = imageNameList.join('/');


      const archived = row.status === 'archived';
      const repoWithArchive = archived ? `${repo}-archive` : repo;

      return html`<a href=https://console.cloud.google.com/artifacts/docker/cpg-common/australia-southeast1/${repoWithArchive}/${image}/sha256:${digest}?project=cpg-common target=_blank>${digest}</a>`;
    }
  }
});

const imageVersionsView = view(imageVersions);
```


<!-- Populate image logs search view -->

```js
// Limit to only getting logs for the first 10 images in the list, otherwise things can get
// a bit out of hand
const versionsMap = imageVersionsView
  .reduce((rr, vv) => {
    rr[vv.image] = rr[vv.image] || [];
    rr[vv.image].push(`'${vv.image}@${vv.digest}'`);
    return rr;
  }, {})

const selectedVersions = Object.values(versionsMap).slice(0, 10).flat().join(',');

const logsForSelectedVersions = await db.query(`
  select
    l.short_path as image,
    i.tags,
    i.status,
    l.timestamp,
    l.principal_email as user,
    l.image_digest as digest
  from logs l
  left join images i
  on i.full_path = l.full_path

  where TRUE
  ${selectedVersions && `and l.image_repository || '/' || l.image_name || '@' || l.image_digest IN (${selectedVersions})`}
  order by timestamp
`)

const logsSearch = Inputs.search(logsForSelectedVersions);
const logsSearchView = view(logsSearch);

```

<!-- Show the log table -->
```js
const logTable = Inputs.table(logsSearchView, {
  layout: 'auto',
  select: false,
  format: {
    timestamp: (dd) => formatTimestamp(dd)
  }
})
```


<!-- Page layout of constructs defined above -->

<div class="grid grid-cols-4">
  <div class="card">
    <h2>Image repositories</h2>
    ${repositoryTable}

  </div>
  <div class="card grid-colspan-3">
    <h2>Images</h2>
    <div style="margin-bottom: 10px;">
      ${imageSearch}
    </div>
    ${imageList}
  </div>

  <div class="card grid-colspan-4">
    <h2>Image versions</h2>
    <div style="margin-bottom: 10px;">
      ${imageVersionSearch}
    </div>
    <div style="margin-bottom: 10px;">
      ${statusSelect}
    </div>
    ${imageVersions}
  </div>

  <div class="card grid-colspan-4">
    <h2>Image version pulls</h2>
    <div style="margin-bottom: 10px;">
      ${logsSearch}
    </div>
    ${logTable}
    <p>
      <em>Logs are only shown for a maximum of 10 selected images</em>
    </p>
  </div>
</div>

<!-- Chart of selected image usage over time -->

<div class="grid grid-cols-1">
  <div class="card">
    ${resize((width) => Plot.plot({
      width: width,
      interactive: true,
      marginLeft: 200,
      x: {
        transform: (d) => new Date(d)
      },
      symbol: {
        legend: true,
        range: ['circle', 'times'],
        domain: ['active', 'archived'],
      },
      marks: [
        Plot.frame(),
        Plot.dot(
          logsSearchView,
          {
            y: 'image',
            x: 'timestamp',
            stroke: 'digest',
            tip: true,
            symbol: 'status',
            channels: {
              tags: 'tags',
              user: 'user',
              status: 'status'
            }
          }
        )
      ]
    }))}
  </div>
</div>
