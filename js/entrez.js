
const https = require('https');
const URL = require('url').URL;
const NAMES = require('../ref_gene_names.json');

class HTTPError extends Error{
  constructor(message,status) {
    super(message);
    this.status = status;
  }
};

class Request {
  constructor(url) {
    this.url = new URL(url);
  }
  async load() {
    const options = {
      hostname: this.url.host,
      port: 443,
      path: this.url.pathname + (this.url.search ? this.url.search : ''),
      method: this.method || 'GET',
      headers: this.headers || {}
    }

    let result = new Promise( (resolve,reject) => {
      let results = [];
      const https_module = https;
      const req = https_module.request(options, (res) => {
        res.on('data', (d) => {
          results.push(d);
        });
        res.on('end', () => {
          if (res.statusCode !== 200) {
            let errmessage = Buffer.concat(results);
            reject(new HTTPError(errmessage,res.statusCode));
            return;
          }
          resolve(Buffer.concat(results))
        });
      });

      req.on('error', reject);

      if (this.body) {
        req.write(this.body);
      }

      req.end();
    });
    return result;
  }
  async loadString() {
    return this.load().then( body => {
      return body.toString();
    });
  }
  async loadJSON() {
    return this.loadString().then( body => {
      return JSON.parse(body);
    });
  }

}

class EntrezMeta {
  init() {
    return Promise.resolve(true);
  }

  lookup(taxid,gene) {
    return NAMES.filter( entry => entry.taxid == taxid && entry.name.toLowerCase() == gene.toLowerCase() )
  }
}

exports = module.exports = new EntrezMeta();